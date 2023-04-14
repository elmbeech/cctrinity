# python libraries
from cc3d.core.PySteppables import *
import numpy as np
import random
import sys


###################
# input variables #
###################

ct = {'WT','wt','wildttype','ROCK1-20','ROCK1','rock1','rock','CDH1','cdh1','cdh'}
ct1 = 'WT'
ct2 = 'WT'
# if ct1 or ct2 is of CDH1 knockdown cell type
ct1_cdh_express_fract = None # {30 .. 90} percentage expression level compared to wt
ct2_cdh_express_fract = None  # {70 .. 90} percentage expression level compared to wt
# perturbation time point
ct1_kd_time = -24  # {-144 .. 120}  # 96
ct2_kd_time = -24 # {-144 .. 120}   # 96
# cell number seeded
ct1_seed = 50
ct2_seed = 50

time_verbose = True
volume_verbose = True
volpsurf_verbose = True
surface_verbose = True
adhesion_verbose = True
velocity_verbose = True

# error check input
if not ct1 in ct:
    sys.exit(f"cell type {ct1} unknown. knowen are {sorted(ct)}.")
if ct1 in {'WT','wt','wildttype'}:
   ct1 = 'wt'
elif ct1 in {'ROCK1-20','ROCK1','rock1','rock'}:
   ct1 = 'ROCK1-20'
elif ct1 in {'CDH1','cdh1','cdh'}:
   if (ct1_cdh_express_fract is None):
       sys.exit(f"Error: ct1 CDH1 expression fraction is {ct1_cdh_express_fract}.")
   ct1 = f'CDH1-{ct1_cdh_express_fract}'

if not ct2 in ct:
    sys.exit(f"cell type {ct2} unknown. knowen are {sorted(ct)}.")
if ct2 in {'WT','wt','wildttype'}:
   ct2 = 'wt'
elif ct2 in {'ROCK1-20','ROCK1','rock1','rock'}:
   ct2 = 'ROCK1-20'
elif ct2 in {'CDH1','cdh1','cdh'}:
   if (ct2_cdh_express_fract is None):
       sys.exit(f"Error: ct2 CDH1 expression fraction is {ct2_cdh_express_fract}.")
   ct2 = f'CDH1-{ct2_cdh_express_fract}'


#################################
# celltype and model parameters #
#################################

# bue 2023-0414: take cell adhesion values form the graner glazer 1992 publication, because they actually work.
# ct-medium the weakest: 16
# ct-ct weak: 14
# ct-ct strong: 2

dd_cellline = {}  # cell mappings
dd_cellline.update({'ROCK1-20': {  # bue: 20 percent rock1 expression level comapred to wild type.
    'adhesion_weak': 8, # graner: ((14 - 2) / 2) + 2 = 9; libby -85  # new data suggests 50% CDH1 expression
    'volume_final': 170,
    #'volume_final_sd': 53.54,
    'volume_edge_final': 228,
    #'volume_edge_final_sd': 60.06,
    #'lambda_volume': 1.0, #  bue: same as init
    'volumepsurface_final': 1.17,  # 1.17,  # was 1.37  # bue: give rock1 kd cells the higher tension because i got rid of dynamic lambda_surface
    'volumepsurface_edge_final': 1.41,  #1.41,  # was 1.97  # bue: give rock 1 kd cells the higher tension because i got rid of dynamic lambda_surface
    #'lambda_surface_final': 1.1,  # wt is 0.5 # bue: set to same as init
    'km_half_time': 46.78,
    'generation_time_final': 20
}})
dd_cellline.update({'wt': {  # bue: wild type expression level
    'adhesion_weak': 2,  # graner: 2; libby -100 # bue: wt has no weak adhesion wt is strong.
    'volume_final': 137,
    #'volume_final_sd': 34.03,  # bue: sd, is this standard deviation?
    'volume_edge_final': 177,
    #'volume_edge_final_sd': 49.66,
    #'lambda_volume': 1.0,
    'volumepsurface_final': 1.12,  # 1.12,
    'volumepsurface_edge_final': 1.32,  # 1.32,
    #'lambda_surface_final': 0.5,  # bue: same as init
    'km_half_time': 48.28, # bue: how can they have km when there is littteraly no protein knocked down? it is just the CDH1 value, that the linear translation works!
    'generation_time_final': 20
}})
dd_cellline.update({'CDH1-0': {   # bue: zero percent cdh1 expression level comapred to wild type
    'adhesion_weak': 14,  # graner: 14; libby -70  # bue: sligthly more than adhesion to medium.
    'volume_final': 123,
    #'volume_final_sd': 51.78,
    'volume_edge_final': 223,
    #'volume_edge_final_sd': 57.96,
    #'lambda_volume': 1.0,
    'volumepsurface_final': 1.10,  # 1.10,
    'volumepsurface_edge_final': 1.23,  # 1.23,
    #'lambda_surface_final': 0.5,  # bue: same as init
    'km_half_time': 48.28,  # time half expression occures after any knockdown.
    'generation_time_final': 18,
}})
# generate mutants for CDH1 (bue: linear correlation assumption)
er_cdh_express_fract = {ct1_cdh_express_fract, ct2_cdh_express_fract}
er_cdh_express_fract.discard(None)
for r_express in er_cdh_express_fract : # bue: beause the same formula is used, sould work for any r_express value between 70 and 90!
    if (r_express < 0.001):
        sys.exit(f'Error: CDH1 expression level {r_express} has to be > 0.001 and should actually be beween in the range 70 .. 90!')
    s_cellline = f'CDH1-{r_express}'
    # initialize cell line
    dd_cellline.update({s_cellline: {
        'adhesion_weak': None,
        'volume_final': None,
        #'volume_final_sd': 51.78,
        'volume_edge_final': None,
        #'volume_edge_final_sd': 57.96,
        #'lambda_volume': 1.0, bue: same as init
        'volumepsurface_final': None,
        'volumepsurface_edge_final': None,
        #'lambda_surface_final': 0.5,  # bue: same as init
        'km_half_time': 48.28,
        'generation_time_final': None,
    }})
    # update initalized cell line
    #r_cdh1 = (100 - r_express) / 100  # fraction expression level compare to wt (relative mutatant decrease)
    r_cdh1 = r_express / 100  # fraction expression level compare to wt
    for s_label in ['adhesion_weak', 'volume_final', 'volume_edge_final', 'volumepsurface_final', 'volumepsurface_edge_final','generation_time_final']:
        #dd_cellline[s_cellline][s_label] = dd_cellline['wt'][s_label] + (dd_cellline['CDH1-100'][s_label] - dd_cellline['wt'][s_label]) * r_cdh1
        dd_cellline[s_cellline][s_label] = dd_cellline['CDH1-0'][s_label] + (dd_cellline['wt'][s_label] - dd_cellline['CDH1-0'][s_label]) * r_cdh1

# do surface calculation
for s_cellline, d_cell in dd_cellline.items():
    d_cell['surface_final'] = d_cell['volume_final'] / d_cell['volumepsurface_final']
    d_cell['surface_edge_final'] = d_cell['volume_edge_final'] / d_cell['volumepsurface_edge_final']

# initialization parameters
volume_init = 157  # volume_edge_init = 177; volume_init = 137
volumepsurface_init = 1.42  # 1.22  # volumepsurface_edge_init = 1.32; volumepsurface_init = 1.12
surface_init = volume_init / volumepsurface_init
adhesion_init = 14  # bue: we don't want that they already start sorting, when they are just kocked down so we set it to the same as weak ct medium and medium medium adhesion
generation_time_init = 20

# constantes
#burn_in_time = 5.0  # to get the seeded system settled in hours?
hill_n = 5.0  # hillpower
lambda_volume = 8.0 # 1.0
lambda_surface = 8.0 # 0.5
lambda_velocity = 9
persistence_velocity = 1

# the morpheus model specifies time in hours
# to map cc3d's mcs to hours we made use of the generation_time and final volumes
# wt only based
# round((137/20 + 177/20) / 2) = 8 #[mcs/h]
# wt, cdh1, and rock1 cell based
# round((137/20 + 177/20 + 170/20 + 228/20 + 123/18 + 223/18) / 6) = 9 #[mcs/h]
# considering the results we settled by 8[mcs/h] or 0.125 [h/mcs]
r_hpmcs = 1/8


#####################
# steppable classes #
#####################

class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)


class GrowthSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # calulate domain center and x60 y20 shape
        i_x = self.dim.x // 2  # floor integer
        i_y = self.dim.y // 2
        i_z = 0
        li_x = list(range(i_x-30, i_x+31))
        li_y = list(range(i_y-10, i_y+11))

        # random seeding
        ct1_drop = ct1_seed
        ct2_drop = ct2_seed
        if (ct1_seed > 0):
            b_flipflop = True
        elif (ct2_seed > 0):
            b_flipflop = False
        else:
            b_flipflop = None
            sys.exit("Error: no ct1 cells ({ct1_seed}) and ct2 cells ({ct2_seed}) set to seed.")
        es_xyz = {None}
        for _ in range(ct1_seed + ct2_seed - 1):
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
            if (b_flipflop):
                self.cell_field[i_x:i_x+2, i_y:i_y+2, i_z] = self.new_cell(self.CT1) # writing to the cell field
                ct1_drop -= 1
                if (ct2_drop > 0):
                    b_flipflop = False
            elif not (b_flipflop):
                self.cell_field[i_x:i_x+2, i_y:i_y+2, i_z] = self.new_cell(self.CT2)  # writing to the cell field
                ct2_drop -= 1
                if (ct1_drop > 0):
                    b_flipflop = True
            else:
                pass

        # set initial constraint.
        for cell in self.cell_list:
            # time
            cell.dict['dt'] = generation_time_init
            cell.dict['dc'] = np.random.uniform(0.0, generation_time_init + r_hpmcs)  # set division time counter
            # volume
            cell.targetVolume = volume_init
            cell.lambdaVolume = lambda_volume
            cell.dict['volumepsurface'] = volumepsurface_init
            # surface
            cell.targetSurface = surface_init
            cell.lambdaSurface = lambda_surface
            # adhesion
            cell.dict['adhesion'] = adhesion_init
            self.adhesionFlexPlugin.setAdhesionMoleculeDensity(cell, 'cdh', cell.dict['adhesion'])
            # velocity
            cell.dict['propulsiveforce'] = lambda_velocity
            cell.dict['alpha'] = persistence_velocity
            cell.dict['angle'] = 2.0 * np.pi * np.random.uniform(-1.0, 1.0)
            cell.dict['dirx'] = np.cos(cell.dict['angle'])
            cell.dict['diry'] = np.sin(cell.dict['angle'])
            cell.lambdaVecX = - cell.dict['propulsiveforce'] * cell.dict['dirx']
            cell.lambdaVecY = - cell.dict['propulsiveforce'] * cell.dict['diry']
            cell.lambdaVecZ = - cell.dict['propulsiveforce'] * cell.dict['diry']

        # plot volumne
        if (volume_verbose):
            self.plot_win_volume = self.add_new_plot_window(
                title='cell volume evolution.',
                x_axis_title='monte_carlo_step',
                y_axis_title='volume',
                x_scale_type='linear',
                y_scale_type='linear',
                grid=True,
                config_options={'legend': True},
            )
            self.plot_win_volume.add_plot("ct1_volume_mean", style='lines', color='cyan', size=2)
            self.plot_win_volume.add_plot("ct2_volume_mean", style='lines', color='blue', size=2)
            self.plot_win_volume.add_plot("ct1_volume_edge_mean", style='lines', color='orange', size=2)
            self.plot_win_volume.add_plot("ct2_volume_edge_mean", style='lines', color='red', size=2)
            self.plot_win_volume.add_plot("ct1_volume", style='dots', color='cyan', size=2)
            self.plot_win_volume.add_plot("ct2_volume", style='dots', color='blue', size=2)
            self.plot_win_volume.add_plot("ct1_volume_edge", style='dots', color='orange', size=2)
            self.plot_win_volume.add_plot("ct2_volume_edge", style='dots', color='red', size=2)
            self.plot_win_volume.add_plot("lambda_volume", style='dots', color='lime', size=2)
            self.plot_win_volume.add_plot("experiment_time", style='dots', color='yellow', size=2)

        # plot surface
        if (surface_verbose):
            self.plot_win_surface = self.add_new_plot_window(
                title='cell surface evolution.',
                x_axis_title='monte_carlo_step',
                y_axis_title='surface',
                x_scale_type='linear',
                y_scale_type='linear',
                grid=True,
                config_options={'legend': True},
            )
            self.plot_win_surface.add_plot("ct1_surface_mean", style='lines', color='cyan', size=2)
            self.plot_win_surface.add_plot("ct2_surface_mean", style='lines', color='blue', size=2)
            self.plot_win_surface.add_plot("ct1_surface_edge_mean", style='lines', color='orange', size=2)
            self.plot_win_surface.add_plot("ct2_surface_edge_mean", style='lines', color='red', size=2)
            self.plot_win_surface.add_plot("ct1_surface", style='dots', color='cyan', size=2)
            self.plot_win_surface.add_plot("ct2_surface", style='dots', color='blue', size=2)
            self.plot_win_surface.add_plot("ct1_surface_edge", style='dots', color='orange', size=2)
            self.plot_win_surface.add_plot("ct2_surface_edge", style='dots', color='red', size=2)
            self.plot_win_surface.add_plot("lambda_surface", style='dots', color='green', size=2)
            self.plot_win_surface.add_plot("experiment_time", style='dots', color='yellow', size=2)

        # plot surface per volume
        if (volpsurf_verbose):
            self.plot_win_volpsurf = self.add_new_plot_window(
                title='cell volume per surface evolution.',
                x_axis_title='monte_carlo_step',
                y_axis_title='volume / surface',
                x_scale_type='linear',
                y_scale_type='linear',
                grid=True,
                config_options={'legend': True},
            )
            self.plot_win_volpsurf.add_plot("ct1_volume_per_surface_mean", style='lines', color='cyan', size=2)
            self.plot_win_volpsurf.add_plot("ct2_volume_per_surface_mean", style='lines', color='blue', size=2)
            self.plot_win_volpsurf.add_plot("ct1_volume_per_surface_edge_mean", style='lines', color='orange', size=2)
            self.plot_win_volpsurf.add_plot("ct2_volume_per_surface_edge_mean", style='lines', color='red', size=2)
            self.plot_win_volpsurf.add_plot("ct1_volumepsurface", style='dots', color='cyan', size=2)
            self.plot_win_volpsurf.add_plot("ct2_volumepsurface", style='dots', color='blue', size=2)
            self.plot_win_volpsurf.add_plot("ct1_volume_per_surface_edge", style='dots', color='orange', size=2)
            self.plot_win_volpsurf.add_plot("ct2_volume_per_surface_edge", style='dots', color='red', size=2)
            self.plot_win_volpsurf.add_plot("experiment_time", style='dots', color='yellow', size=2)

        # plot adhesion
        if (adhesion_verbose):
            self.plot_win_adhesion = self.add_new_plot_window(
                title='cell adhesion evolution.',
                x_axis_title='monte_carlo_step',
                y_axis_title='adhesion',
                x_scale_type='linear',
                y_scale_type='linear',
                grid=True,
                config_options={'legend': True},
            )
            # in
            self.plot_win_adhesion.add_plot("ct1_adhesion", style='dots', color='cyan', size=2)
            self.plot_win_adhesion.add_plot("ct2_adhesion", style='dots', color='blue', size=2)
            self.plot_win_adhesion.add_plot("lambda_adhesion", style='dots', color='lime', size=2)
            # out
            self.plot_win_adhesion.add_plot("total", style='lines', color='maroon', size=2)
            self.plot_win_adhesion.add_plot("medium_medium", style='lines', color='red', size=2)
            self.plot_win_adhesion.add_plot("medium_ct1", style='lines', color='orange', size=2)
            self.plot_win_adhesion.add_plot("medium_ct2", style='lines', color='purple', size=2)
            self.plot_win_adhesion.add_plot("ct1_ct1", style='lines', color='cyan', size=2)
            self.plot_win_adhesion.add_plot("ct1_ct2", style='lines', color='teal', size=2)
            self.plot_win_adhesion.add_plot("ct2_ct2", style='lines', color='blue', size=2)
            # time
            self.plot_win_adhesion.add_plot("experiment_time", style='dots', color='yellow', size=2)

        # plot velocity
        if (velocity_verbose):
            self.plot_win_velocity = self.add_new_plot_window(
                title='cell velocity evolution.',
                x_axis_title='monte_carlo_step',
                y_axis_title='adhesion',
                x_scale_type='linear',
                y_scale_type='linear',
                grid=True,
                config_options={'legend': True},
            )
            self.plot_win_velocity.add_plot("velocity_mean", style='lines', color='purple', size=2)
            self.plot_win_velocity.add_plot("velocity", style='dots', color='magenta', size=2)
            self.plot_win_velocity.add_plot("lambda_velocity", style='dots', color='lime', size=2)
            self.plot_win_velocity.add_plot("persistence_velocity", style='dots', color='cyan', size=2)
            self.plot_win_velocity.add_plot("experiment_time", style='dots', color='yellow', size=2)

        # text output time
        self.msg_win_rudy = self.add_new_message_window(title='Rudy a message to you ...')

    def step(self, mcs):
        # time in hour
        time = (mcs * r_hpmcs) + min(0, ct1_kd_time, ct2_kd_time)
        #print('die zeit:', mcs, time)
        #self.msg_win1.clear()
        self.msg_win_rudy.print(f"simulation time: {mcs}[mcs] {round(time, 3)}[h]")

        # reset variables
        lr_ct1_volume = []
        lr_ct2_volume = []
        lr_ct1_volume_edge = []
        lr_ct2_volume_edge = []
        lr_ct1_surface = []
        lr_ct2_surface = []
        lr_ct1_surface_edge = []
        lr_ct2_surface_edge = []
        lr_ct1_volpsurf = []
        lr_ct2_volpsurf = []
        lr_ct1_volpsurf_edge = []
        lr_ct2_volpsurf_edge = []
        lr_velocity = []

        di_contact = {}
        di_contact.update({'medium_medium': 0})
        di_contact.update({'medium_ct1': 0})
        di_contact.update({'medium_ct2': 0})
        di_contact.update({'ct1_ct1': 0})
        di_contact.update({'ct1_ct2': 0})
        di_contact.update({'ct2_ct2': 0})

        # for each cell
        for cell in self.cell_list:
            # reset variables
            set_volume = None
            set_volumepsurface = None
            set_surface = None
            #set_lambda_volume = None
            #set_lambda_surface = None
            set_adhesion = None

            # get edge or not
            b_edgecell = 0
            for neighbor, _ in self.get_cell_neighbor_data_list(cell):
                if not neighbor:
                    b_edgecell += 1

            #######
            # ct1 #
            #######

            if (cell.type == self.CT1):
                # ct1 adhesion
                if (time < 0):
                    set_adhesion = adhesion_init
                elif (time > ct1_kd_time and cell.dict['adhesion'] >= adhesion_init and cell.dict['adhesion'] <= dd_cellline[ct1]['adhesion_weak']):
                    set_adhesion = (1.0 / (1.0 + (((time - ct1_kd_time) / dd_cellline[ct1]['km_half_time'])**hill_n))) * (adhesion_init - dd_cellline[ct1]['adhesion_weak']) + dd_cellline[ct1]['adhesion_weak']
                elif (time <= ct1_kd_time):
                    set_adhesion = adhesion_init
                else:
                    set_adhesion = dd_cellline[ct1]['adhesion_weak']
                # plot
                if (adhesion_verbose):
                    self.plot_win_adhesion.add_data_point("ct1_adhesion", x=mcs, y=set_adhesion)

                # ct1 lambda_surface
                #if (time < burn_in_time):
                #    set_lambda_surface = lambda_surface_init
                #elif (time > ct1_kd_time and cell.lambdaSurface >= lambda_surface_init and cell.lambdaSurface <= dd_cellline[ct1]['lambda_surface_final']):
                #    set_lambda_surface = (1.0 / (1.0 + (((time - ct1_kd_time) / dd_cellline[ct1]['km_half_time'])**hill_n))) * (lambda_surface_init - dd_cellline[ct1]['lambda_surface_final']) + dd_cellline[ct1]['lambda_surface_final']
                #elif (time <= ct1_kd_time):
                #    set_lambda_surface = lambda_surface_init
                #else:
                #    set_lambda_surface = dd_cellline[ct1]['lambda_surface_final']
                # plot
                #if (surface_verbose):
                #    self.plot_win_surface.add_data_point("ct1_lambda_surface", x=mcs, y=set_lambda_surface)

                # if edge cell
                if (b_edgecell > 0):
                    # ct1_volume_edge
                    if (time < 0):
                        set_volume = volume_init
                    elif (volume_init < dd_cellline[ct1]['volume_edge_final'] and cell.volume < dd_cellline[ct1]['volume_edge_final']):
                        set_volume = volume_init + (dd_cellline[ct1]['volume_edge_final'] - volume_init) * ((time - ct1_kd_time) / 96)
                    else:
                        set_volume = dd_cellline[ct1]['volume_edge_final']
                    # plot
                    if (volume_verbose):
                        self.plot_win_volume.add_data_point("ct1_volume_edge", x=mcs, y=set_volume)
                        lr_ct1_volume_edge.append(cell.volume)

                    # ct1_volume_per_surface_edge and ct1_surface_edge
                    if (time < 0):
                        set_volumepsurface = volumepsurface_init
                    elif (time > ct1_kd_time and cell.dict['volumepsurface'] >= volumepsurface_init and cell.dict['volumepsurface'] <= dd_cellline[ct1]['volumepsurface_edge_final']):
                        set_volumepsurface = (1.0 / (1.0 + (((time - ct1_kd_time) / dd_cellline[ct1]['km_half_time'])**hill_n))) * (volumepsurface_init - dd_cellline[ct1]['volumepsurface_edge_final']) + dd_cellline[ct1]['volumepsurface_edge_final']
                    elif (time <= ct1_kd_time):
                        set_volumepsurface = volumepsurface_init
                    else:
                        set_volumepsurface = dd_cellline[ct1]['volumepsurface_edge_final']
                    set_surface = set_volume / set_volumepsurface
                    # plot
                    if (volpsurf_verbose):
                        self.plot_win_volpsurf.add_data_point("ct1_volume_per_surface_edge", x=mcs, y=set_volumepsurface)
                        lr_ct1_volpsurf_edge.append(set_volumepsurface)
                    if (surface_verbose):
                        self.plot_win_surface.add_data_point("ct1_surface_edge", x=mcs, y=set_surface)
                        lr_ct1_surface_edge.append(cell.surface)

                else:
                    # ct1_volume
                    if (time < 0):
                        set_volume = volume_init
                    elif (volume_init < dd_cellline[ct1]['volume_final'] and cell.volume < dd_cellline[ct1]['volume_final']):
                        set_volume = volume_init + (dd_cellline[ct1]['volume_final'] - volume_init) * ((time - ct1_kd_time) / 96)
                    else:
                        set_volume = dd_cellline[ct1]['volume_final']
                    # plot
                    if (volume_verbose):
                        self.plot_win_volume.add_data_point("ct1_volume", x=mcs, y=set_volume)
                        lr_ct1_volume.append(cell.volume)

                    # ct1_volume_per_surface and ct1_surface
                    if (time < 0):
                        set_volumepsurface = volumepsurface_init
                    elif (time > ct1_kd_time and cell.dict['volumepsurface'] >= volumepsurface_init and cell.dict['volumepsurface'] <= dd_cellline[ct1]['volumepsurface_final']):
                        set_volumepsurface = (1.0 / (1.0 + (((time - ct1_kd_time) / dd_cellline[ct1]['km_half_time'])**hill_n))) * (volumepsurface_init - dd_cellline[ct1]['volumepsurface_final']) + dd_cellline[ct1]['volumepsurface_final']
                    elif (time <= ct1_kd_time):
                        set_volumepsurface = volumepsurface_init
                    else:
                        set_volumepsurface = dd_cellline[ct1]['volumepsurface_final']
                    set_surface = set_volume / set_volumepsurface
                    # plot
                    if (volpsurf_verbose):
                        self.plot_win_volpsurf.add_data_point("ct1_volumepsurface", x=mcs, y=set_volumepsurface)
                        lr_ct1_volpsurf.append(set_volumepsurface)
                    if (surface_verbose):
                        self.plot_win_surface.add_data_point("ct1_surface", x=mcs, y=set_surface)
                        lr_ct1_surface.append(cell.surface)

            #######
            # ct2 #
            #######

            elif (cell.type == self.CT2):
                # ct2_adhesion
                if (time < 0):
                    set_adhesion = adhesion_init
                elif (time > ct2_kd_time and cell.dict['adhesion'] >= adhesion_init and cell.dict['adhesion'] <= dd_cellline[ct2]['adhesion_weak']):
                    set_adhesion = (1.0 / (1.0 + (((time - ct2_kd_time) / dd_cellline[ct2]['km_half_time'])**hill_n))) * (adhesion_init - dd_cellline[ct2]['adhesion_weak']) + dd_cellline[ct2]['adhesion_weak']
                elif (time <= ct2_kd_time):
                    set_adhesion = adhesion_init
                else:
                    set_adhesion = dd_cellline[ct2]['adhesion_weak']
                # plot
                if (adhesion_verbose):
                    self.plot_win_adhesion.add_data_point("ct2_adhesion", x=mcs, y=set_adhesion)

                # ct2_lambda_surface
                #if (time < burn_in_time):
                #    set_lambda_surface = lambda_surface_init
                #elif (time > ct2_kd_time and cell.lambdaSurface >= lambda_surface_init and cell.lambdaSurface <= dd_cellline[ct2]['lambda_surface_final']):
                #    set_lambda_surface = (1.0 / (1.0 + (((time - ct2_kd_time) / dd_cellline[ct1]['km_half_time'])**hill_n))) * (lambda_surface_init - dd_cellline[ct2]['lambda_surface_final']) + dd_cellline[ct2]['lambda_surface_final']
                #elif (time <= ct2_kd_time):
                #    set_lambda_surface = lambda_surface_init
                #else:
                #    set_lambda_surface = dd_cellline[ct2]['lambda_surface_final']
                # plot
                #if (surface_verbose):
                #    self.plot_win_surface.add_data_point("ct2_lambda_surface", x=mcs, y=set_lambda_surface)

                # if edge cell
                if (b_edgecell > 0):
                    # ct2_volume_edge
                    if (time < 0):
                        set_volume = volume_init
                    elif (volume_init < dd_cellline[ct2]['volume_edge_final'] and cell.volume < dd_cellline[ct2]['volume_edge_final']):
                        set_volume = volume_init + (dd_cellline[ct2]['volume_edge_final'] - volume_init) * ((time - ct2_kd_time) / 96)
                    else:
                        set_volume = dd_cellline[ct2]['volume_edge_final']
                    # plot
                    if (volume_verbose):
                        self.plot_win_volume.add_data_point("ct2_volume_edge", x=mcs, y=set_volume)
                        lr_ct2_volume_edge.append(cell.volume)

                    # ct2_volume_per_surface_edge and ct2_surface_edge
                    if (time < 0):
                        set_volumepsurface = volumepsurface_init
                    elif (time > ct2_kd_time and cell.dict['volumepsurface'] >= volumepsurface_init and cell.dict['volumepsurface'] <= dd_cellline[ct2]['volumepsurface_edge_final']):
                        set_volumepsurface = (1.0 / (1.0 + (((time - ct2_kd_time) / dd_cellline[ct2]['km_half_time'])**hill_n))) * (volumepsurface_init - dd_cellline[ct2]['volumepsurface_edge_final']) + dd_cellline[ct2]['volumepsurface_edge_final']
                    elif (time <= ct2_kd_time):
                        set_volumepsurface = volumepsurface_init
                    else:
                        set_volumepsurface = dd_cellline[ct2]['volumepsurface_edge_final']
                    set_surface = set_volume / set_volumepsurface
                    # plot
                    if (volpsurf_verbose):
                        self.plot_win_volpsurf.add_data_point("ct2_volume_per_surface_edge", x=mcs, y=set_volumepsurface)
                        lr_ct2_volpsurf_edge.append(set_volumepsurface)
                    if (surface_verbose):
                        self.plot_win_surface.add_data_point("ct2_surface_edge", x=mcs, y=set_surface)
                        lr_ct2_surface_edge.append(cell.surface)

                else:
                    # ct2_volume
                    if (time < 0):
                        set_volume = volume_init
                    elif (volume_init < dd_cellline[ct2]['volume_final'] and cell.volume < dd_cellline[ct2]['volume_final']):
                        set_volume = volume_init + (dd_cellline[ct2]['volume_final'] - volume_init) * ((time - ct2_kd_time) / 96)
                    else:
                        set_volume = dd_cellline[ct2]['volume_final']
                    # plot
                    if (volume_verbose):
                        self.plot_win_volume.add_data_point("ct2_volume", x=mcs, y=set_volume)
                        lr_ct2_volume.append(cell.volume)

                    # ct2_volume per surface
                    if (time < 0):
                        set_volumepsurface = volumepsurface_init
                    elif (time > ct2_kd_time and cell.dict['volumepsurface'] >= volumepsurface_init and cell.dict['volumepsurface'] <= dd_cellline[ct2]['volumepsurface_final']):
                        set_volumepsurface = (1.0 / (1.0 + (((time - ct2_kd_time) / dd_cellline[ct2]['km_half_time'])**hill_n))) * (volumepsurface_init - dd_cellline[ct2]['volumepsurface_final']) + dd_cellline[ct2]['volumepsurface_final']
                    elif (time <= ct2_kd_time):
                        set_volumepsurface = volumepsurface_init
                    else:
                        set_volumepsurface = dd_cellline[ct2]['volumepsurface_final']
                    set_surface = set_volume / set_volumepsurface
                    # plot
                    if (volpsurf_verbose):
                        self.plot_win_volpsurf.add_data_point("ct2_volumepsurface", x=mcs, y=set_volumepsurface)
                        lr_ct2_volpsurf.append(set_volumepsurface)
                    if (surface_verbose):
                        self.plot_win_surface.add_data_point("ct2_surface", x=mcs, y=set_surface)
                        lr_ct2_surface.append(cell.surface)


            ############
            # ct error #
            ############

            else:
                sys.exit("unknowen cell type detetced: {cell.type}")

            ##########################
            # update cell parameters #
            ##########################

            # volume
            cell.targetVolume = int(set_volume)
            #cell.lambdaVolume = lambda_volume_int  # set_lambda_volume  bue: lambda_volume undergoes no calcualtion.
            # surface
            cell.dict['volumepsurface'] = set_volumepsurface
            cell.targetSurface = int(set_surface)
            #cell.lambdaSurface = set_lambda_surface
            # adhesion
            cell.dict['adhesion'] = set_adhesion
            self.adhesionFlexPlugin.setAdhesionMoleculeDensity(cell, 'cdh', cell.dict['adhesion'])
            # velocity
            cell.dict['angle'] = cell.dict['angle'] + 2.0 * np.pi * np.random.uniform(-1.0,1.0) * (1 - cell.dict['alpha'])
            cell.dict['dirx'] = np.cos(cell.dict['angle'])
            cell.dict['diry'] = np.sin(cell.dict['angle'])
            cell.lambdaVecX = - cell.dict['propulsiveforce'] * cell.dict['dirx']
            cell.lambdaVecY = - cell.dict['propulsiveforce'] * cell.dict['diry']
            cell.lambdaVecZ = 0.0

            #####################
            # calculate contact #
            #####################
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                i_contact = common_surface_area / 2
                if neighbor:
                    #print(f"cell.id: {cell.id}\tcell.type: {cell.type}\tneighbor.id: {neighbor.id}\tneighbor.type: {neighbor.type}\tcommon_surface_area: {common_surface_area}")
                    if (cell.type == 1) and (neighbor.type == 1):
                        di_contact['ct1_ct1'] += i_contact

                    elif ((cell.type == 1) and (neighbor.type == 2)) or ((cell.type == 2) and (neighbor.type == 1)):
                        di_contact['ct1_ct2'] += i_contact

                    elif (cell.type == 2) and (neighbor.type == 2):
                        di_contact['ct2_ct2'] += i_contact

                    else:
                        sys.exit(f"Error: unknowen cell.type neighbor.type detected: {cell.type } {neighbor.type}")
                else:
                    #print(f"cell.type: {cell.type}\tcell.type: {cell.type}\tmedium common_surface_area: {common_surface_area}\t{round(i_contact, 3)}")
                    if (cell.type == 1):
                        di_contact['medium_ct1'] += common_surface_area

                    elif (cell.type == 2):
                        di_contact['medium_ct2'] += common_surface_area

                    else:
                        sys.exit(f"Error: unknowen cell.id with medium contact detected: {cell.type}")
            #print(f"cell.id: {cell.id}\t{sorted(di_contact.items())}")

        # plot means constantes and time
        # volume
        if (volume_verbose):
            self.plot_win_volume.add_data_point("ct1_volume_mean", x=mcs, y=np.mean(lr_ct1_volume))
            self.plot_win_volume.add_data_point("ct2_volume_mean", x=mcs, y=np.mean(lr_ct2_volume))
            self.plot_win_volume.add_data_point("ct1_volume_edge_mean", x=mcs, y=np.mean(lr_ct1_volume_edge))
            self.plot_win_volume.add_data_point("ct2_volume_edge_mean", x=mcs, y=np.mean(lr_ct2_volume_edge))
            self.plot_win_volume.add_data_point("lambda_volume", x=mcs, y=lambda_volume)
            self.plot_win_volume.add_data_point("experiment_time", x=mcs, y=time)
        # surface
        if (surface_verbose):
            self.plot_win_surface.add_data_point("ct1_surface_mean", x=mcs, y=np.mean(lr_ct1_surface))
            self.plot_win_surface.add_data_point("ct2_surface_mean", x=mcs, y=np.mean(lr_ct2_surface))
            self.plot_win_surface.add_data_point("ct1_surface_edge_mean", x=mcs, y=np.mean(lr_ct1_surface_edge))
            self.plot_win_surface.add_data_point("ct2_surface_edge_mean", x=mcs, y=np.mean(lr_ct2_surface_edge))
            self.plot_win_surface.add_data_point("lambda_surface", x=mcs, y=lambda_surface)
            self.plot_win_surface.add_data_point("experiment_time", x=mcs, y=time)
        if (volpsurf_verbose):
            self.plot_win_volpsurf.add_data_point("ct1_volume_per_surface_mean", x=mcs, y=np.mean(lr_ct1_volpsurf))
            self.plot_win_volpsurf.add_data_point("ct2_volume_per_surface_mean", x=mcs, y=np.mean(lr_ct2_volpsurf))
            self.plot_win_volpsurf.add_data_point("ct1_volume_per_surface_edge_mean", x=mcs, y=np.mean(lr_ct1_volpsurf_edge))
            self.plot_win_volpsurf.add_data_point("ct2_volume_per_surface_edge_mean", x=mcs, y=np.mean(lr_ct2_volpsurf_edge))
            self.plot_win_volpsurf.add_data_point("experiment_time", x=mcs, y=time)
        # adhesion
        if (adhesion_verbose):
            # in
            lambda_adhesion = self.get_xml_element("lamdaAdhesion").cdata
            self.plot_win_adhesion.add_data_point("lambda_adhesion", x=mcs, y=lambda_adhesion)
            # out
            self.plot_win_adhesion.add_data_point("total", x=mcs, y=sum(di_contact.values()))
            self.plot_win_adhesion.add_data_point("medium_ct1", x=mcs, y=di_contact['medium_ct1'] )
            self.plot_win_adhesion.add_data_point("medium_ct2", x=mcs, y=di_contact['medium_ct2'])
            self.plot_win_adhesion.add_data_point("ct1_ct1", x=mcs, y=di_contact['ct1_ct1'])
            self.plot_win_adhesion.add_data_point("ct1_ct2", x=mcs, y=di_contact['ct1_ct2'])
            self.plot_win_adhesion.add_data_point("ct2_ct2", x=mcs, y=int(di_contact['ct2_ct2']))
            # time
            self.plot_win_adhesion.add_data_point("experiment_time", x=mcs, y=time)
            #print(f"{mcs}[mcs]\t{time}[h]\ttotal: {sum(di_contact.values())}\t{sorted(di_contact.items())}")
        # velocity
        if (velocity_verbose):
            self.plot_win_velocity.add_data_point("lambda_velocity", x=mcs, y=lambda_velocity)
            self.plot_win_velocity.add_data_point("persistence_velocity", x=mcs, y=persistence_velocity)
            self.plot_win_velocity.add_data_point("experiment_time", x=mcs, y=time)

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


    def start(self):
        # plot window time
        if (time_verbose):
            self.plot_win_time = self.add_new_plot_window(
                title='time evolution',
                x_axis_title='monte_carlo_step',
                y_axis_title='hour',
                x_scale_type='linear',
                y_scale_type='linear',
                grid=True,
                config_options={'legend': True},
            )
            self.plot_win_time.add_plot("experiment_time", style='dots', color='yellow', size=2)
            self.plot_win_time.add_plot("generation_time", style='dots', color='red', size=2)
            self.plot_win_time.add_plot("division_time", style='dots', color='orange', size=2)


    def step(self, mcs):

        ########
        # time #
        ########

        # time in hour
        time = (mcs * r_hpmcs) + min(0, ct1_kd_time, ct2_kd_time)

        # get cells to devide
        cells_to_divide=[]
        for cell in self.cell_list:
            set_generation_time = None


            #######
            # ct1 #
            #######

            if (cell.type == self.CT1):
                # ct1 gernation time
                if (time < 0):
                    set_generation_time = generation_time_init
                elif (time > ct1_kd_time):
                    set_generation_time = (1.0 / (1.0 + (((time - ct1_kd_time) / dd_cellline[ct1]['km_half_time'])**hill_n))) * (generation_time_init - dd_cellline[ct1]['generation_time_final']) + dd_cellline[ct1]['generation_time_final']
                elif (time <= ct1_kd_time):
                    set_generation_time = generation_time_init
                else:
                    set_generation_time = dd_cellline[ct1]['generation_time_final']


            #######
            # ct2 #
            #######

            elif (cell.type == self.CT2):
                # ct2_generation_time
                if (time < 0):
                    set_generation_time = generation_time_init
                elif (time > ct2_kd_time):
                    set_generation_time = (1.0 / (1.0 + (((time - ct2_kd_time) / dd_cellline[ct2]['km_half_time'])**hill_n))) * (generation_time_init - dd_cellline[ct2]['generation_time_final']) + dd_cellline[ct2]['generation_time_final']
                elif (time <= ct2_kd_time):
                    set_generation_time = generation_time_init
                else:
                    set_generation_time = dd_cellline[ct2]['generation_time_final']


            ############
            # ct error #
            ############

            else:
                 sys.exit("unknowen cell type detetced: {cell.type}")


            ##########################
            # update cell parameters #
            ##########################

            cell.dict['dt'] = set_generation_time
            if (cell.dict['dc'] >= cell.dict['dt']):
                cells_to_divide.append(cell)
            if (time >= 0):
                cell.dict['dc'] += r_hpmcs

            # plot cell division time
            if (time_verbose) and (mcs % 8):
                self.plot_win_time.add_data_point("division_time", x=mcs, y=cell.dict['dc'])

        # plot experiment time in hour and cell generation time
        if (time_verbose) and (mcs % 8):
            self.plot_win_time.add_data_point("experiment_time", x=mcs, y=time)
            self.plot_win_time.add_data_point("generation_time", x=mcs, y=cell.dict['dt'])


        #################
        # cell division #
        #################

        for cell in cells_to_divide:
            cell.dict['dt'] = generation_time_init
            cell.dict['dc'] = 0  # set divison time counter to zero
            self.divide_cell_random_orientation(cell)

    def update_attributes(self):
        self.parent_cell.targetVolume /= 2.0  # reducing parent target volume
        self.clone_parent_2_child()

