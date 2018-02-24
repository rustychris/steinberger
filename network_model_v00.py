import numpy as np

from stompy.spatial import wkb2shp, field
from shapely import geometry
from stompy.grid import unstructured_grid
from stompy.plot import plot_wkb
from stompy import utils

from stompy import memoize

##
dem=field.GdalGrid('steinberger-dem-crop.tif')

@memoize.memoize()
def poly_to_depths(poly):
    depths=self.dem.F[ self.dem.polygon_mask(poly) ]
    depths.sort()
    return depths

## 
class RoutingModel(object):
    eta_range=[-2,2]
    N_discrete_eta=100

    T=0 # simulation time in seconds
    # Set up the mass conservation problem
    dt=600 # seconds

    def __init__(self,bounds_shp,flow_lines_shp,dem):
        self.bounds_shp=bounds_shp
        self.flow_lines_shp=flow_lines_shp

        self.load_dem(dem)
        self.load_shps()
        self.init_geometry()
        self.instrument_grid()
        self.init_flow_lines()
        self.extract_ordering()
        self.preprocess_hypsometry()
        self.preprocess_barycentric()

        self.init_particles()

    def load_dem(self,dem):
        if isinstance(dem,str):
            dem=field.GdalGrid(dem)

        self.dem=dem
    def load_shps(self):
        self.bounds=wkb2shp.shp2geom(self.bounds_shp)
        self.flow_lines=wkb2shp.shp2geom(self.flow_lines_shp)

    def init_geometry(self):
        """
        Compile a grid with simple polygon cells from the 
        boundaries in bounds_shp
        """
        g=unstructured_grid.UnstructuredGrid()

        for geo in self.bounds['geom']:
            pnts=np.array(geo)

            A=pnts[:-1]
            B=pnts[1:]

            for a,b in zip(A,B):
                na=g.add_or_find_node(x=a)
                nb=g.add_or_find_node(x=b)
                try:
                    j=g.add_edge(nodes=[na,nb])
                except g.GridException:
                    pass

        cycles=g.find_cycles(max_cycle_len=g.Nnodes())

        polys=[geometry.Polygon( g.nodes['x'][cycle] )
               for cycle in cycles ]

        gc=unstructured_grid.UnstructuredGrid(max_sides=max( [len(cyc) for cyc in cycles] ))

        gc.copy_from_grid(g)

        for cyc in cycles:
            gc.add_cell(nodes=cyc)

        self.g=gc

    def instrument_grid(self):
        # 0: not a flow edge
        # 1: 'downstream' flow is c1 to c2
        # -1: 'downstream' flow is c2 to c1
        self.g.add_edge_field('sign',np.zeros(self.g.Nedges(),'i4'),on_exists='overwrite')
        self.g.edges['sign']=0 #

        self.fluxes=np.zeros(g.Nedges(), 'f8')

    def init_flow_lines(self):
        # lateral boundary conditions
        self.edge_sources=[] # [ (name, edge_index), ... ]
        self.edge_sinks=[]   # [ (name, edge_index), ... ]

        # flow_lines are used to add downstream cell relationships
        for i in range(len(self.flow_lines)):
            line=self.flow_lines['geom'][i]
            name=self.flow_lines['name'][i]

            # Quick and dirty, find the order in which this line visits cells,
            # include entering/leaving the grid
            dx=10.
            samples=[np.array(line.interpolate(d)) for d in np.arange(0,line.length,dx)]
            cells=np.array( [max(-1,self.g.select_cells_nearest(s,inside=True))
                             for s in samples] )

            for samp_i in range(1,len(samples)):
                a=cells[samp_i-1]
                b=cells[samp_i]
                if a==b:
                    continue

                # would like to identify the specific edge and mark it.
                samp_line=geometry.LineString( [ samples[samp_i-1], samples[samp_i]] )
                hits=self.g.select_edges_intersecting(samp_line)
                j=np.nonzero(hits)[0][0]

                if self.g.edges['cells'][j,0]==a:
                    self.g.edges['sign'][j]=1
                else:
                    self.g.edges['sign'][j]=-1

                if a<0:
                    self.edge_sources.append( [name,j] )
                if b<0:
                    self.edge_sinks.append( [name,j])


    def downstream(self,j):
        s=self.g.edges['sign'][j]
        if s==0:
            return -1 # not a flow edge
        return self.g.edges['cells'][j,(1+s)//2]

    def upstream(self,j):
        s=self.g.edges['sign'][j]
        if s==0:
            return -1 # not a flow edge
        return self.g.edges['cells'][j,(1-s)//2]

    def downstream_edges(self,c):
        js=[]
        for j in self.g.cell_to_edges(c):
            if self.upstream(j)==c:
                js.append(j)
        return js

    def upstream_edges(self,c):
        js=[]
        for j in self.g.cell_to_edges(c):
            if self.downstream(j)==c:
                js.append(j)
        return js

    def extract_ordering(self):
        # get the upstream to downstream ordering:
        visited=np.zeros(self.g.Ncells(),'i4') - 1
        
        def visit_edge(j,count):
            c=self.upstream(j)
            if c>=0:
                return visit_cell(c,count)
            else:
                return count

        def visit_cell(c,count):
            if visited[c]>=0:
                return count
            else:
                visited[c]=count
                count+=1
                for j in self.upstream_edges(c):
                    count=visit_edge(j,count)
                return count

        count=0    
        for name,j in edge_sinks:
            count=visit_edge(j,count)

        self.visited=visited
        self.cell_order=np.argsort( -visited )[ :np.sum(visited>=0) ]

    def preprocess_hypsometry(self):
        """
        Discretize the depth distributions to make later volume requests
        fast.
        """
        discrete_eta=np.linspace(self.eta_range[0],self.eta_range[1],self.N_discrete_eta)
        self.discrete_eta=discrete_eta
        self.preprocess_hypsometry_cells()
        self.preprocess_hypsometry_edges()

    def preprocess_hypsometry_cells(self):
        # want a function that efficiently maps an eta to a volume 
        # Evaluate this for each eta in discrete_eta, after which
        # it will be quick to interpolate
        vol_for_eta=np.zeros( (self.g.Ncells(),len(discrete_eta)), 'f8')
        Apixel=self.dem.dx * self.dem.dy

        for c in range(self.g.Ncells()):
            print("Hypsometry for cell %d"%c)
            poly=self.g.cell_polygon(c)
            depths=poly_to_depths(poly)
            cutoffs=np.searchsorted(depths,discrete_eta)
            for eta_i,eta in enumerate(discrete_eta):
                if cutoffs[eta_i]==0:
                    vol_for_eta[c,eta_i]=0.0
                else:
                    vol_for_eta[c,eta_i]=np.sum( (eta-depths[:cutoffs[eta_i]]) * Apixel )
        self.vol_for_eta=vol_for_eta

    def preprocess_hypsometry_edges(self):
        # Edges
        fluxarea_for_eta=np.zeros( (self.g.Nedges(), len(self.discrete_eta)), 'f8')
        dx=1.0
        e2c=self.g.edge_to_cells()
        e_len=self.g.edges_length()

        for j in range(self.g.Nedges()):

            if self.g.edges['sign'][j]==0:
                continue
            # print("Depth for j=%d"%j)

            c1,c2=e2c[j]
            if (c1>=0) and (c2>=0):
                j_a=(e2c[:,0]==c1) & (e2c[:,1]==c2)
                j_b=(e2c[:,0]==c2) & (e2c[:,1]==c1)
                all_j=np.nonzero( j_a|j_b )[0]
            else:
                all_j=[j]

            all_depths=[]

            print("  %d sub-edges"%len(all_j))

            for jj in all_j:
                pnts=self.g.nodes['x'][self.g.edges['nodes'][jj]]
                dist=e_len[jj]
                n_samps=int(round(dist/dx))
                alpha=np.linspace(0,1,n_samps+1)
                alpha=0.5*(alpha[:-1]+alpha[1:])
                alpha=alpha[:,None] # for broadcasting
                samps=(1-alpha)*pnts[0] + alpha*pnts[1]
                depths=self.dem(samps)
                all_depths.append(depths)
            depths=np.concatenate(all_depths)
            depths.sort()
            cutoffs=np.searchsorted(depths,self.discrete_eta)
            for eta_i,eta in enumerate(self.discrete_eta):
                if cutoffs[eta_i]==0:
                    fluxarea_for_eta[j,eta_i]=0.0
                else:
                    fluxarea_for_eta[j,eta_i]=np.sum( (eta-depths[:cutoffs[eta_i]]) * dx )
            self.fluxarea_for_eta=fluxarea_for_eta

    def preprocess_barycentric(self):
        # get a nice mapping of cells to active edges:
        active_edges=self.g.edges['sign']!=0

        cell_to_active_edges=[]
        cell_to_active_signs=[]

        for c in range(self.g.Ncells()):
            # edges for which a negative flux is *INTO* this cell
            neg_edges=np.nonzero( (self.g.edges['cells'][:,0]==c) & active_edges )[0]
            # edges for which a positive flux is *INTO* this cell
            pos_edges=np.nonzero( (self.g.edges['cells'][:,1]==c) & active_edges )[0]
            cell_to_active_edges.append( np.concatenate([neg_edges,pos_edges]) )
            cell_to_active_signs.append( np.concatenate([-np.ones_like(neg_edges),
                                                         np.ones_like(pos_edges)]) )

        # max_active_edges=max( [len(x) for x in cell_to_active_edges] )
        self.cell_to_active_edges=cell_to_active_edges
        self.cell_to_active_signs=cell_to_active_signs

    def init_particles(self):
        self.inactive_particles=[]
        self.particles=[] # list of active particle locations [ (cell,[coords]), ... ]

    def add_particle(self,loc):
        self.particles.append( loc )

    def cell_eta(self,t):
        omega=2*np.pi/(12*3600)
        return np.cos(omega*t) * np.ones(self.g.Ncells(),'f8')
    def edge_eta(self,t):
        cell_eta=self.cell_eta(t)
        
        c1=self.g.edges['cells'][:,0].copy()
        c2=self.g.edges['cells'][:,1].copy()
        c2[c2<0]=c1[c2<0]
        c1[c1<0]=c2[c1<0]
        # choose the higher eta:
        return np.maximum(cell_eta[c1],cell_eta[c2])

    def volumes(self,t):
        etas=self.cell_eta(t)

        assert etas.min()>=self.discrete_eta[0]
        assert etas.max()<=self.discrete_eta[-1]

        vols=np.zeros(self.g.Ncells(),'f8')
        for c in range(self.g.Ncells()):

            vols[c] = np.interp(etas[c],
                                self.discrete_eta,self.vol_for_eta[c,:])
        return vols
    def fluxareas(self,t):
        etas=self.edge_eta(t)

        assert etas.min()>=self.discrete_eta[0]
        assert etas.max()<=self.discrete_eta[-1]

        areas=np.zeros(self.g.Nedges(),'f8')
        for j in range(self.g.Nedges()):
            areas[j] = np.interp(etas[j],
                                 self.discrete_eta,self.fluxarea_for_eta[j,:])
        return areas

    def step(self):
        """
        advance from self.T to self.T+self.dt 
        """
        self.step_flows()
        self.step_particles()
        self.T+=self.dt

    def step_flows(self):
        # starting volume:
        self.vol0=vol0=self.volumes(t=self.T)
        # volume at end of step
        self.vol1=vol1=self.volumes(t=self.T+self.dt)

        # rates, m3/s, positive is left-to-right across edge
        self.fluxes=fluxes=np.zeros(self.g.Nedges(),'f8') 

        # Source fluxes are set as dirichlet BCs on flow rate
        # m3/s -- should come from inputs, positive is downstream/source
        flow_rate=10.0

        for source_name,j in self.edge_sources:
            fluxes[j]+=self.g.edges['sign'][j]*flow_rate

        # Starting with the most upstream cells, route fluxes downstream
        for c in self.cell_order:
            v=vol0[c]
            for j in self.upstream_edges(c):
                v+=self.dt*(fluxes[j]*self.g.edges['sign'][j])
            
            surplus_v=v - vol1[c]
            surplus_flux=surplus_v / self.dt

            js=self.downstream_edges(c)
            assert len(js)==1,"only drainage trees are supported at this point"
            for j in js:
                fluxes[j]+=surplus_flux * self.g.edges['sign'][j]
                
    # Particle methods:
    # Particle location are described by a cell index, and barycentric coordinates
    def particle_input(self,j_input):
        """
        place a particle on this edge, input to the downstream cell
        """
        c_input=self.downstream(j_input)
        coord=[int(j_input==j)
               for j in self.cell_to_active_edges[c_input]]
        return [c_input,np.array(coord)]

    def bary_to_xy(self,loc):
        c,coord=loc

        centers=self.g.edges_center()[self.cell_to_active_edges[c]]
        return (centers*coord[:,None]).sum(axis=0)

    def step_particles(self):
        i=0
        while i<len(self.particles):
            loc=self.particles[i]
            new_loc=self.particle_advect(loc,self.dt)
            if new_loc[0]>=0: # still active
                self.particles[i]=new_loc
                i+=1
            else:
                self.inactive_particles.append(loc)
                del self.particles[i]

    def particle_advect(self,loc,dt):
        c,coord=loc
        if c<0:
            # out of domain
            return loc
        js=self.cell_to_active_edges[c]
        flux_in=self.fluxes[js]*self.cell_to_active_signs[c]
        flux_adv=flux_in-flux_in.mean() # remove convergence 
        V=0.5*(self.vol0[c]+self.vol1[c])
        norm_delta=flux_adv/V

        t_step=dt
        j_cross=None
        for i in range(len(coord)):
            if norm_delta[i]>=0:
                continue # edge has in-flux, moving away from it
            t_to_cross=(1.0-coord[i])/(-norm_delta[i])
            if t_to_cross<t_step:
                t_step=t_to_cross
                j_cross=js[i]

        if j_cross is not None:
            c1,c2=self.g.edges['cells'][j_cross]
            if c1==c:
                new_c=c2
            elif c2==c:
                new_c=c1
            else:
                assert False,"how did that happen?"
            new_coord=[ int(new_j)==j_cross
                        for new_j in self.cell_to_active_edges[new_c] ]

            return self.particle_advect((new_c,new_coord), dt-t_step)
        else:
            new_coord=coord - t_step*norm_delta

        assert np.abs( new_coord.sum()-1 ) < 1e-5
        assert new_coord.min()>=0.0
        assert new_coord.max()<=1.0

        return (c,new_coord)
    
bounds_shp='finite_volume-bounds.shp'
flow_lines_shp="flow_lines_v01.shp"

router=RoutingModel(bounds_shp,flow_lines_shp,dem=dem)

for src in router.edge_sources:
    loc=router.particle_input( src[1] )
    router.add_particle(loc)

## 

router.step()

self=router

edge_mask=self.g.edges['sign']!=0

plt.figure(1).clf()
fig,ax=plt.subplots(num=1)

self.g.plot_cells(color='0.7',lw=0.5,mask=self.visited<0,ax=ax)
ccoll=self.g.plot_cells(values=self.visited,lw=0.5,mask=self.visited>=0,ax=ax,alpha=0.3)
ccoll.set_edgecolors('k')
ccoll.set_cmap('jet')

[plot_wkb.plot_wkb(geo,color='k',ax=ax) for geo in self.flow_lines['geom']]

fluxes=self.fluxes
areas=self.fluxareas(t=0).clip(1,np.inf)
self.g.plot_edges( labeler=lambda j,r: "%.1f (%.2fm/s)"%(self.g.edges['sign'][j]*fluxes[j],
                                                         self.g.edges['sign'][j]*fluxes[j]/areas[j]),
                   mask=edge_mask, ax=ax )

p_xy=np.zeros( (len(self.particles),2), np.float64)

for p_i,loc in enumerate(self.particles):
    p_xy[p_i,:]=router.bary_to_xy(loc) 

ax.plot(p_xy[:,0],p_xy[:,1],'bo')

## 

loc=router.particle_input(j_input)
xy=router.bary_to_xy(loc)

plt.plot([xy[0]],[xy[1]],'ro')

pnts=[]
for t in np.arange(0,500*3600,10000):
    new_loc=router.particle_advect(self,loc,t)
    if new_loc[0]<0:
        break
    new_xy=router.bary_to_xy(new_loc)
    pnts.append(new_xy)
pnts=np.array(pnts)

plt.plot(pnts[:,0],pnts[:,1],'g-o')

## 

# not quite right with the barycentric advection.
# currently written as the weighted sum of the nodes
# so (1,0,0) means exactly at node 0.
# this is breaking when a particle enters at, say, node 0,
# and there is flow in from node 1, and out of 2.

# the setup: 
#  N flows
#  V(t) volume
#  

# would this work better with reaches and junctions?
# then location would be linear within a reach, and
# probabilistic at a junction proportional to flow rates.
# that would be more representative of the actual configuration.
# as it stands, a reach with two upstream and one downstream
# links will have some funny behavior -- traveling from one
# upstream to the other requires exchanging the entire volume
# of the reach.
