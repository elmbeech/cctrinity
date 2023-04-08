# library
from cc3d.core.PySteppables import *
import numpy as np
import random
import sys

# suggest user input variable
ct = {'wt','wildttype','ROCK1-20','ROCK1','rock1','rock','CDH1','cdh1','cdh'}
ct1 = 'wt'
ct2 = 'ROCK1'
# if ct1 or ct2 is CDH1 knockdown celline
ct1_cdh_kd_fract = None  # {70 .. 90}
ct2_cdh_kd_fract = None  # {70 .. 90}
# perturbation time point
ct1_kd_time = -144  # {-144 .. 96}
ct2_kd_time = -74  # {-144 .. 96}
# cell number seeded
ct1_seed = 50
ct2_seed = 50

# error check and update input
if not ct1 in ct:
    sys.exit(f"cell type {ct1} unknown. knowen are {sorted(ct)}.")
if ct1 in {'wt','wildttype'}:
   ct1 = 'wt'
elif ct1 in {'ROCK1-20','ROCK1','rock1','rock'}:
   ct1 = 'ROCK1-20'
elif ct1 in {'CDH1','cdh1','cdh'}:
   ct1 = f'CDH1-{ct1_cdh_kd_fract}'

if not ct2 in ct:
    sys.exit(f"cell type {ct2} unknown. knowen are {sorted(ct)}.")
if ct2 in {'wt','wildttype'}:
   ct2 = 'wt'
elif ct2 in {'ROCK1-20','ROCK1','rock1','rock'}:
   ct2 = 'ROCK1-20'
elif ct2 in {'CDH1','cdh1','cdh'}:
   ct2 = f'CDH1-{ct2_cdh_kd_fract}'

# <Variable symbol="ct1_adhesion" value="-100"/>
# <Variable symbol="ct2_adhesion" value="-100"/>
# <Variable symbol="ct1_ct2" value="max(ct1_adhesion,ct2_adhesion)"/>
# <Variable symbol="ctall_medium_adhesion" value="0"/>

# const and var
dd_cellline = {}  # cell mappings
dd_cellline.update({'wt': {
    'adhesion_weak': -100,  # bue: why is this negative?
    'volume_final': 137,
    #'volume_final_sd': 34.03,  # bue: sd, is this standard deviation?
    'volume_edge_final': 177,
    #'volume_edge_final_sd': 49.66,
    #'lambda_volume': 1.0,
    'volumepsurface_final': 1.12,
    'volumepsurface_edge_final': 1.32,
    'lambda_surface_final': 0.5,
    'km_half_time': 48.28, # bue: how can they have km when there is littteraly no protein knocked down?
    'generation_time_final': 20
}})
dd_cellline.update({'ROCK1-20': {
    'adhesion_weak': -85, # new data suggests 50% CDH1 expression
    'volume_final': 170,
    #'volume_final_sd': 53.54,
    'volume_edge_final': 228,
    #'volume_edge_final_sd': 60.06,
    #'lambda_volume': 1.0,
    'volumepsurface_final': 1.17,  # was 1.37
    'volumepsurface_edge_final': 1.41,  # was 1.97
    'lambda_surface_final': 1.1,  # wt is 0.5
    'km_half_time': 46.78,
    'generation_time_final': 20
}})
dd_cellline.update({'CDH1-0': {
    'adhesion_weak': -70,
    'volume_final': 123,
    #'volume_final_sd': 51.78,
    'volume_edge_final': 223,
    #'volume_edge_final_sd': 57.96,
    #'lambda_volume': 1.0,
    'volumepsurface_final': 1.10,
    'volumepsurface_edge_final': 1.23,
    'lambda_surface_final': 0.5,
    'km_half_time': 48.28,  # time half expression occures after any knockdown.
    'generation_time_final': 18,
}})
# generate mutants for CDH1
er_cdh_kd_fract = {ct1_cdh_kd_fract, ct2_cdh_kd_fract}
er_cdh_kd_fract.discard(None)
for r_kd in er_cdh_kd_fract : # percent of cdh1 gene knowdown # bue: beause the same formula is used, sould work for any r_kd value between 70 and 90!
    if (r_kd < 0.001):
        sys.exit(f'Error: CDH1 knockdown percentage {r_kd} has to be > 0.001 and should actually be beween in the range 70 .. 90!')
    s_cellline = f'CDH1-{r_kd}'
    # initialize cell line
    dd_cellline.update({s_cellline: {
        'adhesion_weak': None,
        'volume_final': None,
        #'volume_final_sd': 51.78,
        'volume_edge_final': None,
        #'volume_edge_final_sd': 57.96,
        #'lambda_volume': 1.0,
        'volumepsurface_final': None,
        'volumepsurface_edge_final': None,
        'lambda_surface_final': 0.5,
        'km_half_time': 48.28,
        'generation_time_final': 20,
    }})
    # update initalized cell line
    r_cdh1 = (100 - r_kd) / 100  # relative mutatant decrease
    for s_label in ['adhesion_weak', 'volume_final', 'volume_edge_final', 'volumepsurface_final', 'volumepsurface_edge_final']:
        dd_cellline[s_cellline][s_label] = dd_cellline['wt'][s_label] + (dd_cellline['CDH1-0'][s_label] - dd_cellline['wt'][s_label]) * r_cdh1

# do surface calculation
for s_cellline, d_cell in dd_cellline.items():
    d_cell['surface_final'] = d_cell['volume_final'] / d_cell['volumepsurface_final']
    d_cell['surface_edge_final'] = d_cell['volume_edge_final'] / d_cell['volumepsurface_edge_final']


# initialize parameter
volume_init = 157  # volume_edge_init = 177; volume_init = 137
lambda_volume_int = 1.0
volumepsurface_init = 1.22  # volumepsurface_edge_init = 1.32; volumepsurface_init = 1.12
surface_init = volume_init / volumepsurface_init
lambda_surface_init = 0.5
adhesion_init = -100  # bue: we don't want that they already start sorting, when they are just kocked down.
persistent_motion_str_init  = 9
generation_time_init = 20

# constantes
burn_in_time = 5.0  # to get the seeded system settled in hours
hill_n = 5.0

# the morpheus model specifies time in hours
# to map cc3d's mcs to hours we made use of the generation_time and final volumes
# wt only based
# round((137/20 + 177/20) / 2) = 8 #[mcs/h]
# wt, cdh1, and rock1 cell based
# round((137/20 + 177/20 + 170/20 + 228/20 + 123/18 + 223/18) / 6) = 9 #[mcs/h]
# considering the results we settled by 8[mcs/h] or 0.125 [h/mcs]
r_hpmcs = 1/8


# steppable classes
class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)


class GrowthSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        ## seeding cells
        # calx domain center calx
        i_x = self.dim.x // 2  # floor integer
        i_y = self.dim.y // 2
        i_z = 0
        li_x = list(range(i_x-30, i_x+31))
        li_y = list(range(i_y-10, i_y+11))

        # random seed in domain x60 y20
        es_xyz = {None}
        for _ in range(ct1_seed):
            s_xyz0 = None
            s_xyz1 = None
            s_xyz2 = None
            s_xyz3 = None
            while (s_xyz0 in es_xyz) or (s_xyz1 in es_xyz) or (s_xyz2 in es_xyz) or (s_xyz3 in es_xyz):
                i_x = random.choice(li_x)
                i_y = random.choice(li_y)
                s_xyz0 = f'{i_x}{i_y}{i_z}'
                s_xyz2 = f'{i_x}{i_y+1}{i_z}'
                s_xyz1 = f'{i_x+1}{i_y}{i_z}'
                s_xyz3 = f'{i_x+1}{i_y+1}{i_z}'
            es_xyz = es_xyz.union({s_xyz0, s_xyz1, s_xyz2, s_xyz3})
            self.cell_field[i_x:i_x+2, i_y:i_y+2, i_z] = self.new_cell(self.CT1) # writing to the cell field

        for _ in range(ct2_seed):
            s_xyz0 = None
            s_xyz1 = None
            s_xyz2 = None
            s_xyz3 = None
            while (s_xyz0 in es_xyz) or (s_xyz1 in es_xyz) or (s_xyz2 in es_xyz) or (s_xyz3 in es_xyz):
                i_x = random.choice(li_x)
                i_y = random.choice(li_y)
                s_xyz0 = f'{i_x}{i_y}{i_z}'
                s_xyz2 = f'{i_x}{i_y+1}{i_z}'
                s_xyz1 = f'{i_x+1}{i_y}{i_z}'
                s_xyz3 = f'{i_x+1}{i_y+1}{i_z}'
            es_xyz = es_xyz.union({s_xyz0, s_xyz1, s_xyz2, s_xyz3})
            self.cell_field[i_x:i_x+2, i_y:i_y+2, i_z] = self.new_cell(self.CT2)  # writing to the cell field

        # set initial volume and surface constraints
        for cell in self.cell_list:
            cell.dict['dt'] = generation_time_init
            cell.dict['dc'] = 0  # set division time cunter to zero
            cell.targetVolume = volume_init
            cell.lambdaVolume = lambda_volume_int
            cell.dict['volumepsurface'] = volumepsurface_init
            cell.targetSurface = surface_init
            cell.lambdaSurface = lambda_surface_init
            cell.dict['adhesion'] = adhesion_init

    def step(self, mcs):
        # get time in hour
        time = (mcs / r_hpmcs) + min(0, ct1_kd_time, ct2_kd_time)

        # for each cell
        for cell in self.cell_list:
            set_ct1_adhesion = None  # i have to set all combinations!?
            set_ct2_adhesion = None
            set_volume = None
            set_volumepsurface = None
            set_surface = None
            set_lambda_volume = None
            set_lambda_surface = None
            set_adhesion = None

            # get edge cell or not
            isEdgeCell = 0
            for neighbor, _ in self.get_cell_neighbor_data_list(cell):
                if not neighbor:
                    isEdgeCell += 1

            ### ct1 ###
            if (cell.type == self.CT1):
                # ct1 adhesion
                if (time < burn_in_time):
                    set_adhesion = adhesion_init
                elif (time > ct1_kd_time and cell.dict['adhesion'] >= adhesion_init and cell.dict['adhesion'] <= dd_cellline[ct1]['adhesion_weak']):
                    set_adhesion = (1.0 / (1.0 + (((time - ct1_kd_time) / dd_cellline[ct1]['km_half_time'])**hill_n))) * (adhesion_init - dd_cellline[ct1]['adhesion_weak']) + dd_cellline[ct1]['adhesion_weak']
                elif (time <= ct1_kd_time):
                    set_adhesion = adhesion_init
                else:
                    set_adhesion = dd_cellline[ct1]['adhesion_weak']

                # ct1 lambda_surface
                if (time < burn_in_time):
                    set_lambda_surface = lambda_surface_init
                elif (time > ct1_kd_time and cell.lambdaSurface >= lambda_surface_init and cell.lambdaSurface <= dd_cellline[ct1]['lambda_surface_final']):
                    set_lambda_surface = (1.0 / (1.0 + (((time - ct1_kd_time) / dd_cellline[ct1]['km_half_time'])**hill_n))) * (lambda_surface_init - dd_cellline[ct1]['lambda_surface_final']) + dd_cellline[ct1]['lambda_surface_final']
                elif (time <= ct1_kd_time):
                    set_lambda_surface = lambda_surface_init
                else:
                    set_lambda_surface = dd_cellline[ct1]['lambda_surface_final']

                # if edge cell
                if (isEdgeCell > 0):
                    # ct volume_edge
                    if (time < burn_in_time):
                        set_volume = volume_init
                    elif (volume_init < dd_cellline[ct1]['volume_edge_final'] and cell.volume < dd_cellline[ct1]['volume_edge_final']):
                        set_volume = volume_init + (dd_cellline[ct1]['volume_edge_final'] - volume_init) * ((time - ct1_kd_time) / 96)
                    else:
                        set_volume = dd_cellline[ct1]['volume_edge_final']

                    #ct1_volume per surface_edge
                    if (time < burn_in_time):
                        set_volumepsurface = volumepsurface_init
                    elif (time > ct1_kd_time and cell.dict['volumepsurface'] >= volumepsurface_init and cell.dict['volumepsurface'] <= dd_cellline[ct1]['volumepsurface_edge_final']):
                        set_volumepsurface = (1.0 / (1.0 + (((time - ct1_kd_time) / dd_cellline[ct1]['km_half_time'])**hill_n))) * (volumepsurface_init - dd_cellline[ct1]['volumepsurface_edge_final']) + dd_cellline[ct1]['volumepsurface_edge_final']
                    elif (time <= ct1_kd_time):
                        set_volumepsurface = volumepsurface_init
                    else:
                        set_volumepsurface = dd_cellline[ct1]['volumepsurface_edge_final']

                else:
                    # ct1 volume
                    if (time < burn_in_time):
                        set_volume = volume_init
                    elif (volume_init < dd_cellline[ct1]['volume_final'] and cell.volume < dd_cellline[ct1]['volume_final']):
                        set_volume = volume_init + (dd_cellline[ct1]['volume_final'] - volume_init) * ((time - ct1_kd_time) / 96)
                    else:
                        set_volume = dd_cellline[ct1]['volume_final']

                    # ct1 volume per surface
                    if (time < burn_in_time):
                        set_volumepsurface = volumepsurface_init
                    elif (time > ct1_kd_time and cell.dict['volumepsurface'] >= volumepsurface_init and cell.dict['volumepsurface'] <= dd_cellline[ct1]['volumepsurface_final']):
                        set_volumepsurface = (1.0 / (1.0 + (((time - ct1_kd_time) / dd_cellline[ct1]['km_half_time'])**hill_n))) * (volumepsurface_init - dd_cellline[ct1]['volumepsurface_final']) + dd_cellline[ct1]['volumepsurface_final']
                    elif (time <= ct1_kd_time):
                        set_volumepsurface = volumepsurface_init
                    else:
                        set_volumepsurface = dd_cellline[ct1]['volumepsurface_final']

            ### ct2 ###
            elif (cell.type == self.CT2):
                # ct2_adhesion
                if (time < burn_in_time):
                    set_adhesion = adhesion_init
                elif (time > ct2_kd_time and cell.dict['adhesion'] >= adhesion_init and cell.dict['adhesion'] <= dd_cellline[ct2]['adhesion_weak']):
                    set_adhesion = (1.0 / (1.0 + (((time - ct2_kd_time) / dd_cellline[ct2]['km_half_time'])**hill_n))) * (adhesion_init - dd_cellline[ct2]['adhesion_weak']) + dd_cellline[ct2]['adhesion_weak']
                elif (time <= ct2_kd_time):
                    set_ct2_adhesion = adhesion_init
                else:
                    set_adhesion = dd_cellline[ct2]['adhesion_weak']

                # ct2_lambda_surface
                if (time < burn_in_time):
                    set_lambda_surface = lambda_surface_init
                elif (time > ct2_kd_time and cell.lambdaSurface >= lambda_surface_init and cell.lambdaSurface <= dd_cellline[ct2]['lambda_surface_final']):
                    set_lambda_surface = (1.0 / (1.0 + (((time - ct2_kd_time) / dd_cellline[ct1]['km_half_time'])**hill_n))) * (lambda_surface_init - dd_cellline[ct2]['lambda_surface_final']) + dd_cellline[ct2]['lambda_surface_final']
                elif (time <= ct2_kd_time):
                    set_lambda_surface = lambda_surface_init
                else:
                    set_lambda_surface = dd_cellline[ct2]['lambda_surface_final']

                # if edge cell
                if (isEdgeCell > 0):
                    # ct2_volume_edge
                    if (time < burn_in_time):
                        set_volume = volume_init
                    elif (volume_init < dd_cellline[ct2]['volume_edge_final'] and cell.volume < dd_cellline[ct2]['volume_edge_final']):
                        set_volume = volume_init + (dd_cellline[ct2]['volume_edge_final'] - volume_init) * ((time - ct2_kd_time) / 96)
                    else:
                        set_volume = dd_cellline[ct2]['volume_edge_final']

                    # ct2_volume per surface_edge
                    if (time < burn_in_time):
                        set_volumepsurface = volumepsurface_init
                    elif (time > ct2_kd_time and cell.dict['volumepsurface'] >= volumepsurface_init and cell.dict['volumepsurface'] <= dd_cellline[ct2]['volumepsurface_edge_final']):
                        set_volumepsurface = (1.0 / (1.0 + (((time - ct2_kd_time) / dd_cellline[ct2]['km_half_time'])**hill_n))) * (volumepsurface_init - dd_cellline[ct2]['volumepsurface_edge_final']) + dd_cellline[ct2]['volumepsurface_edge_final']
                    elif (time <= ct2_kd_time):
                        set_volumepsurface = volumepsurface_init
                    else:
                        set_volumepsurface = dd_cellline[ct2]['volumepsurface_edge_final']

                else:
                    # ct2_volume
                    if (time < burn_in_time):
                        set_volume = volume_init
                    elif (volume_init < dd_cellline[ct2]['volume_final'] and cell.volume < dd_cellline[ct2]['volume_final']):
                        set_volume = volume_init + (dd_cellline[ct2]['volume_final'] - volume_init) * ((time - ct2_kd_time) / 96)
                    else:
                        set_volume = dd_cellline[ct2]['volume_final']

                    # ct2_volume per surface
                    if (time < burn_in_time):
                        set_volumepsurface = volumepsurface_init
                    elif (time > ct2_kd_time and cell.dict['volumepsurface'] >= volumepsurface_init and cell.dict['volumepsurface'] <= dd_cellline[ct2]['volumepsurface_final']):
                        set_volumepsurface = (1.0 / (1.0 + (((time - ct2_kd_time) / dd_cellline[ct2]['km_half_time'])**hill_n))) * (volumepsurface_init - dd_cellline[ct2]['volumepsurface_final']) + dd_cellline[ct2]['volumepsurface_final']
                    elif (time <= ct2_kd_time):
                        set_volumepsurface = volumepsurface_init
                    else:
                        set_volumepsurface = dd_cellline[ct2]['volumepsurface_final']

            ### ct error ###
            else:
                sys.exit("unknowen cell type detetced: {cell.type}")

            ### update cell ###
            cell.targetVolume = int(set_volume)
            cell.lambdaVolume = lambda_volume_int  # set_lambda_volume  bue: lambda_volume undergoes no calcualtion.
            cell.dict['volumepsurface'] = set_volumepsurface
            cell.targetSurface = int(set_volume / set_volumepsurface)
            cell.lambdaSurface = set_lambda_surface
            cell.dict['adhesion'] = set_adhesion

            # bue: todo
            # set all cell-cell adession values, so that they have effect.

            # bue: todo persistent motion
            # <PropertyVector symbol="dv" value="0.0, 0.0, 0.0" name="direction_vector"/>
            # <PersistentMotion protrusion="true" decay-time="0.003" strength="9"/>


    def finish(self):
        """
        Called after the last MCS to wrap up the simulation. Good place to close files and do pos    t-processing
        """
        pass


    def on_stop(self):
        """
        Called if the simulation is stopped before the last MCS
        """
        self.finish()


class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)

    def step(self, mcs):
        # get time in hour
        time = (mcs / r_hpmcs) + min(0, ct1_kd_time, ct2_kd_time)

        # get cells to devide
        cells_to_divide=[]
        for cell in self.cell_list:
            set_generation_time = None

            ### ct1 ###
            if (cell.type == self.CT1):
                # ct1 gernation time
                if (time < burn_in_time):
                    set_generation_time = generation_time_init
                elif (time > ct1_kd_time):
                    set_generation_time = (1.0 / (1.0 + (((time - ct1_kd_time) / dd_cellline[ct1]['km_half_time'])**hill_n))) * (generation_time_init - dd_cellline[ct1]['generation_time_final']) + dd_cellline[ct1]['generation_time_final']
                elif (time <= ct1_kd_time):
                    set_generation_time = generation_time_init
                else:
                    set_generation_time = dd_cellline[ct1]['generation_time_final']

            ### ct2 ###
            elif (cell.type == self.CT2):
                # ct2_generation_time
                if (time < burn_in_time):
                    set_generation_time = generation_time_init
                elif (time > ct2_kd_time):
                    set_generation_time = (1.0 / (1.0 + (((time - ct2_kd_time) / dd_cellline[ct2]['km_half_time'])**hill_n))) * (generation_time_init - dd_cellline[ct2]['generation_time_final']) + dd_cellline[ct2]['generation_time_final']
                elif (time <= ct2_kd_time):
                    set_generation_time = generation_time_init
                else:
                    set_generation_time = dd_cellline[ct2]['generation_time_final']

            ### ct error ###
            else:
                 sys.exit("unknowen cell type detetced: {cell.type}")

            ### update cell ###
            cell.dict['dt'] = set_generation_time
            if (cell.dict['dc'] >= cell.dict['dt']):
                cells_to_divide.append(cell)
            if (time >= 0):
                cell.dict['dc'] += r_hpmcs

        for cell in cells_to_divide:
            cell.dict['dt'] = generation_time_init
            cell.dict['dc'] = 0  # set divison time counter to zero
            self.divide_cell_random_orientation(cell)

    def update_attributes(self):
        self.parent_cell.targetVolume /= 2.0  # reducing parent target volume
        self.clone_parent_2_child()

