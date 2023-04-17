# python libraries
from cc3d.core.PySteppables import *
import numpy as np
import random
import sys


#################################
# celltype and model parameters #
#################################

# bue 2023-0414: take cell adhesion values form the graner glazer 1992 publication, because they actually work.
# ct-medium the weakest: 16
# ct-ct weak: 14
# ct-ct strong: 2
adhesion_max = 2

dd_cellline = {}  # cell mappings
dd_cellline.update({'WT': {  # bue: wild type expression level
    'adhesion_min': adhesion_max,  # graner: 2; libby -100 # bue: wt has no weak adhesion wt is strong.
    'volume_final': 137,
    'volume_edge_final': 177,
    'volumepsurface_final': 1.12,
    'volumepsurface_edge_final': 1.32,
    'km_half_time': 48.28, # bue: how can they have km when there is littteraly no protein knocked down? it is just the CDH1 value, that the linear translation works!
    'generation_time_final': 20
}})
dd_cellline.update({'CDH1-0': {   # bue: zero percent cdh1 expression level comapred to wild type
    'adhesion_min': 28,  # graner: 14; libby -70  # bue: sligthly more than adhesion to medium.
    'volume_final': 123,
    'volume_edge_final': 223,
    'volumepsurface_final': 1.10,
    'volumepsurface_edge_final': 1.23,
    'km_half_time': 48.28,  # time half expression occures after any knockdown.
    'generation_time_final': 18,
}})
dd_cellline.update({'ROCK1-20': {  # bue: 20 percent rock1 expression level comapred to wild type.
    'adhesion_min': 16, # graner: ((14 - 2) / 2) + 2 = 8; libby -85  # new data suggests 50% CDH1 expression
    'volume_final': 170,
    'volume_edge_final': 228,
    'volumepsurface_final': 1.17,  # was 1.37  # bue: give rock1 kd cells the higher tension because i got rid of dynamic lambda_surface
    'volumepsurface_edge_final': 1.41,  # was 1.97  # bue: give rock 1 kd cells the higher tension because i got rid of dynamic lambda_surface
    'km_half_time': 46.78,
    'generation_time_final': 20
}})

# constantes
adhesion_init = 32  # bue: we don't want that they already start sorting, when they are just kocked down so we set it to the same as weak ct medium and medium medium adhesion
hill_n = 5.0  # hillpower
lambda_volume = 4.0 # 1.0
lambda_surface = 4.0 # 0.5
lambda_velocity = 9  # 9
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

class GrowthSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        pass

    def add_steering_panel(self):
        # define input variables, specify cell type
        self.add_steering_param(name='ct1', val='WT', enum=['WT','CDH1','ROCK1-20'], widget_name='combobox')
        self.add_steering_param(name='ct2', val='WT', enum=['WT','CDH1','ROCK1-20'], widget_name='combobox')
        # if ct1 or ct2 is of CDH1 knockdown cell type {70 .. 90} percentage expression level compared to wt
        self.add_steering_param(name='ct1_cdh_express_fract', val=100, min_val=0, max_val=100, widget_name='slider')
        self.add_steering_param(name='ct2_cdh_express_fract', val=100, min_val=0, max_val=100, widget_name='slider')
        # knock down time {-144 .. 120}  # 96
        self.add_steering_param(name='ct1_kd_time', val=0, min_val=-96, max_val=96, widget_name='slider')
        self.add_steering_param(name='ct2_kd_time', val=0, min_val=-96, max_val=96, widget_name='slider')
        # seeding
        self.add_steering_param(name='ct1_seed', val=50, min_val=0, max_val=100, widget_name='slider')
        self.add_steering_param(name='ct2_seed', val=50, min_val=0, max_val=100, widget_name='slider')
        # plotting windows
        self.add_steering_param(name='time_plot', val='off', enum=['off','on'], widget_name='pull-down')
        self.add_steering_param(name='volume_plot', val='off', enum=['off','on'], widget_name='pull-down')
        self.add_steering_param(name='volpsurf_plot', val='off', enum=['off','on'], widget_name='pull-down')
        self.add_steering_param(name='surface_plot', val='off', enum=['off','on'], widget_name='pull-down')
        self.add_steering_param(name='adhesion_plot', val='on', enum=['off','on'], widget_name='pull-down')
        self.add_steering_param(name='contact_plot', val='on', enum=['off','on'], widget_name='pull-down')
        self.add_steering_param(name='velocity_plot', val='off', enum=['off','on'], widget_name='pull-down')
        # should I stay or should I go?
        self.add_steering_param(name='processing', val='WAIT!', enum=['WAIT!','GO!'], widget_name='pull-down')

    def process_steering_panel_data(self):
        # cell type
        global ct1
        global ct2
        ct1 = self.get_steering_param('ct1')
        ct2 = self.get_steering_param('ct2')
        # if ct1 or ct2 is of CDH1 cell type percentage expression level compared to wt
        ct1_cdh_express_fract = self.get_steering_param('ct1_cdh_express_fract')
        ct2_cdh_express_fract = self.get_steering_param('ct2_cdh_express_fract')
        if ct1 == 'CDH1':
            ct1 = f'CDH1-{ct1_cdh_express_fract}'
        if ct2 == 'CDH1':
            ct2 = f'CDH1-{ct2_cdh_express_fract}'
        # knockdown time
        global ct1_kd_time
        global ct2_kd_time
        ct1_kd_time = self.get_steering_param('ct1_kd_time')
        ct2_kd_time = self.get_steering_param('ct2_kd_time')
        # plotwindows
        global time_verbose
        global volume_verbose
        global volpsurf_verbose
        global surface_verbose
        global adhesion_verbose
        global contact_verbose
        global velocity_verbose
        time_verbose = False
        volume_verbose = False
        volpsurf_verbose = False
        surface_verbose = False
        adhesion_verbose = False
        contact_verbose = False
        velocity_verbose = False
        if self.get_steering_param('time_plot') == 'on':
            time_verbose = True
        if self.get_steering_param('volume_plot') == 'on':
            volume_verbose = True
        if self.get_steering_param('volpsurf_plot') == 'on':
            volpsurf_verbose = True
        if self.get_steering_param('surface_plot') == 'on':
            surface_verbose = True
        if self.get_steering_param('adhesion_plot') == 'on':
            adhesion_verbose = True
        if self.get_steering_param('contact_plot') == 'on':
            contact_verbose = True
        if self.get_steering_param('velocity_plot') == 'on':
            velocity_verbose = True

        # generate mutants for CDH1 (bue: linear correlation assumption)
        er_cdh_express_fract = {ct1_cdh_express_fract, ct2_cdh_express_fract}
        er_cdh_express_fract.discard(None)
        for r_express in er_cdh_express_fract : # bue: beause the same formula is used, sould work for any r_express value between 70 and 90!
            if (r_express > 0.001):
                s_cellline = f'CDH1-{r_express}'
                # initialize cell line
                dd_cellline.update({s_cellline: {
                    'adhesion_min': None,
                    'volume_final': None,
                    'volume_edge_final': None,
                    'volumepsurface_final': None,
                    'volumepsurface_edge_final': None,
                    'km_half_time': 48.28,  # same as CDH1-0
                    'generation_time_final': None,
                }})
                # update initalized cell line
                r_cdh1 = r_express / 100  # fraction expression level compare to wt
                for s_label in ['adhesion_min', 'volume_final', 'volume_edge_final', 'volumepsurface_final', 'volumepsurface_edge_final','generation_time_final']:
                    dd_cellline[s_cellline][s_label] = dd_cellline['CDH1-0'][s_label] + (dd_cellline['WT'][s_label] - dd_cellline['CDH1-0'][s_label]) * r_cdh1

        # do surface calculation
        for s_cellline, d_cell in dd_cellline.items():
            d_cell['surface_final'] = d_cell['volume_final'] / d_cell['volumepsurface_final']
            d_cell['surface_edge_final'] = d_cell['volume_edge_final'] / d_cell['volumepsurface_edge_final']

    def step(self, mcs):
        if (mcs==0):
            # hold
            while self.get_steering_param('processing') != 'GO!':
               pass
            self.process_steering_panel_data()
            ct1_seed = self.get_steering_param('ct1_seed')
            ct2_seed = self.get_steering_param('ct2_seed')

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
                # ct1
                if (b_flipflop):
                    self.cell_field[i_x:i_x+2, i_y:i_y+2, i_z] = self.new_cell(self.CT1) # writing to the cell field
                    # next
                    ct1_drop -= 1
                    if (ct2_drop > 0):
                        b_flipflop = False
                # ct2
                elif not (b_flipflop):
                    self.cell_field[i_x:i_x+2, i_y:i_y+2, i_z] = self.new_cell(self.CT2)  # writing to the cell field
                    # next
                    ct2_drop -= 1
                    if (ct1_drop > 0):
                        b_flipflop = True
                # pass
                else:
                    pass
    
            # initialize cell parameters
            for cell in self.cell_list:
                # ct1
                if (cell.type == self.CT1):
                    # time
                    cell.dict['dt'] = dd_cellline[ct1]['generation_time_final']
                    cell.dict['dc'] = np.random.uniform(0.0, dd_cellline[ct1]['generation_time_final'] + r_hpmcs)  # set division time counter
                    # volume
                    cell.targetVolume = dd_cellline[ct1]['volume_final']
                    cell.lambdaVolume = lambda_volume
                    # surface
                    cell.targetSurface = dd_cellline[ct1]['volume_final'] / dd_cellline[ct1]['volumepsurface_final']
                    cell.lambdaSurface = lambda_surface
                # ct2
                elif (cell.type == self.CT2):
                    # time
                    cell.dict['dt'] = dd_cellline[ct2]['generation_time_final']
                    cell.dict['dc'] = np.random.uniform(0.0, dd_cellline[ct2]['generation_time_final'] + r_hpmcs)  # set division time counter
                    # volume
                    cell.targetVolume = dd_cellline[ct2]['volume_final']
                    cell.lambdaVolume = lambda_volume
                    # surface
                    cell.targetSurface = dd_cellline[ct2]['volume_final'] / dd_cellline[ct2]['volumepsurface_final']
                    cell.lambdaSurface = lambda_surface
                # ct error
                else:
                    sys.exit(f"unknowen cell type detetced: {cell.type}")
                # adhesion
                self.adhesionFlexPlugin.setAdhesionMoleculeDensity(cell, 'cdh', adhesion_init)
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
                self.plot_win_adhesion.add_plot("medium_adhesion", style='dots', color='maroon', size=2)
                self.plot_win_adhesion.add_plot("ct1_adhesion", style='dots', color='cyan', size=2)
                self.plot_win_adhesion.add_plot("ct2_adhesion", style='dots', color='blue', size=2)
                self.plot_win_adhesion.add_plot("wt_adhesion", style='dots', color='teal', size=2)
                self.plot_win_adhesion.add_plot("lambda_adhesion", style='dots', color='lime', size=2)
                # time
                self.plot_win_adhesion.add_plot("experiment_time", style='dots', color='yellow', size=2)
    
            # plot contact
            if (contact_verbose):
                self.plot_win_contact = self.add_new_plot_window(
                    title='cell contact evolution.',
                    x_axis_title='monte_carlo_step',
                    y_axis_title='contact length',
                    x_scale_type='linear',
                    y_scale_type='linear',
                    grid=True,
                    config_options={'legend': True},
                )
                self.plot_win_contact.add_plot("total", style='lines', color='maroon', size=2)
                self.plot_win_contact.add_plot("medium_medium", style='lines', color='red', size=2)
                self.plot_win_contact.add_plot("medium_ct1", style='lines', color='orange', size=2)
                self.plot_win_contact.add_plot("medium_ct2", style='lines', color='purple', size=2)
                self.plot_win_contact.add_plot("ct1_ct1", style='lines', color='cyan', size=2)
                self.plot_win_contact.add_plot("ct1_ct2", style='lines', color='teal', size=2)
                self.plot_win_contact.add_plot("ct2_ct2", style='lines', color='blue', size=2)
                # time
                self.plot_win_contact.add_plot("experiment_time", style='dots', color='yellow', size=2)
    
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
            #self.msg_win_rudy = self.add_new_message_window(title='Rudy a message to you ...')  # bue: yet unavilable on nanohub.

        # time in hour
        time = (mcs * r_hpmcs) + min(0, ct1_kd_time, ct2_kd_time)
        #print('die zeit:', mcs, time)
        #self.msg_win1.clear()  # bue: yet unavilable on nanohub.
        #self.msg_win_rudy.print(f"simulation time: {mcs}[mcs] {round(time, 3)}[h]")  # bue: yet unavilable on nanohub.

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
        #lr_ct1_adhesion_mean = []
        #lr_ct2_adhesion_mean = []
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
                else:
                    set_adhesion = dd_cellline[ct1]['adhesion_min'] + ((adhesion_max - dd_cellline[ct1]['adhesion_min']) * (time - ct1_kd_time)**hill_n) / (dd_cellline[ct1]['km_half_time']**hill_n + (time - ct1_kd_time)**hill_n)
                # plot
                if (adhesion_verbose):
                    self.plot_win_adhesion.add_data_point("ct1_adhesion", x=mcs, y=set_adhesion)

                # if edge cell
                if (b_edgecell > 0):
                    # ct1_volume_edge
                    if (time < 0):
                        set_volume = dd_cellline[ct1]['volume_edge_final']
                    else:
                        set_volume = dd_cellline[ct1]['volume_edge_final'] + ((dd_cellline['WT']['volume_edge_final'] - dd_cellline[ct1]['volume_edge_final']) * (time - ct1_kd_time)**hill_n) / (dd_cellline[ct1]['km_half_time']**hill_n + (time - ct1_kd_time)**hill_n)
                    # plot
                    if (volume_verbose):
                        self.plot_win_volume.add_data_point("ct1_volume_edge", x=mcs, y=set_volume)
                        lr_ct1_volume_edge.append(cell.volume)

                    # ct1_volume_per_surface_edge and ct1_surface_edge
                    if (time < 0):
                        set_volumepsurface = dd_cellline[ct1]['volumepsurface_edge_final']
                    else:
                        set_volumepsurface = dd_cellline[ct1]['volumepsurface_edge_final'] + ((dd_cellline['WT']['volumepsurface_edge_final'] - dd_cellline[ct1]['volumepsurface_edge_final']) * (time - ct1_kd_time)**hill_n) / (dd_cellline[ct1]['km_half_time']**hill_n + (time - ct1_kd_time)**hill_n)
                    # calculate surface
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
                        set_volume = dd_cellline[ct1]['volume_final']
                    else:
                        set_volume = dd_cellline[ct1]['volume_final'] + ((dd_cellline['WT']['volume_final'] - dd_cellline[ct1]['volume_final']) * (time - ct1_kd_time)**hill_n) / (dd_cellline[ct1]['km_half_time']**hill_n + (time - ct1_kd_time)**hill_n)
                    # plot
                    if (volume_verbose):
                        self.plot_win_volume.add_data_point("ct1_volume", x=mcs, y=set_volume)
                        lr_ct1_volume.append(cell.volume)

                    # ct1_volume_per_surface and ct1_surface
                    if (time < 0):
                        set_volumepsurface = dd_cellline[ct1]['volumepsurface_final']
                    else:
                        set_volumepsurface = dd_cellline[ct1]['volumepsurface_final'] + ((dd_cellline['WT']['volumepsurface_final'] - dd_cellline[ct1]['volumepsurface_final']) * (time - ct1_kd_time)**hill_n) / (dd_cellline[ct1]['km_half_time']**hill_n + (time - ct1_kd_time)**hill_n)
                    # calculate surface
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
                else:
                    set_adhesion = dd_cellline[ct2]['adhesion_min'] + ((adhesion_max - dd_cellline[ct2]['adhesion_min']) * (time - ct2_kd_time)**hill_n) / (dd_cellline[ct2]['km_half_time']**hill_n + (time - ct2_kd_time)**hill_n)
                # plot
                if (adhesion_verbose):
                    self.plot_win_adhesion.add_data_point("ct2_adhesion", x=mcs, y=set_adhesion)

                # if edge cell
                if (b_edgecell > 0):
                    # ct2_volume_edge
                    if (time < 0):
                        set_volume = dd_cellline[ct2]['volume_edge_final']
                    else:
                        set_volume = dd_cellline[ct2]['volume_edge_final'] + ((dd_cellline['WT']['volume_edge_final'] - dd_cellline[ct2]['volume_edge_final']) * (time - ct2_kd_time)**hill_n) / (dd_cellline[ct2]['km_half_time']**hill_n + (time - ct2_kd_time)**hill_n)
                    # plot
                    if (volume_verbose):
                        self.plot_win_volume.add_data_point("ct2_volume_edge", x=mcs, y=set_volume)
                        lr_ct2_volume_edge.append(cell.volume)

                    # ct2_volume_per_surface_edge and ct2_surface_edge
                    if (time < 0):
                        set_volumepsurface = dd_cellline[ct2]['volumepsurface_edge_final']
                    else:
                        set_volumepsurface = dd_cellline[ct2]['volumepsurface_edge_final'] + ((dd_cellline['WT']['volumepsurface_edge_final'] - dd_cellline[ct2]['volumepsurface_edge_final']) * (time - ct2_kd_time)**hill_n) / (dd_cellline[ct2]['km_half_time']**hill_n + (time - ct2_kd_time)**hill_n)
                    # calculate surface
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
                        set_volume = dd_cellline[ct2]['volume_final']
                    else:
                        set_volume = dd_cellline[ct2]['volume_final'] + ((dd_cellline['WT']['volume_final'] - dd_cellline[ct2]['volume_final']) * (time - ct2_kd_time)**hill_n) / (dd_cellline[ct2]['km_half_time']**hill_n + (time - ct2_kd_time)**hill_n)
                    # plot
                    if (volume_verbose):
                        self.plot_win_volume.add_data_point("ct2_volume", x=mcs, y=set_volume)
                        lr_ct2_volume.append(cell.volume)

                    # ct2_volume per surface
                    if (time < 0): #or (time <= ct2_kd_time):
                        set_volumepsurface = dd_cellline[ct2]['volumepsurface_final']
                    else:
                        set_volumepsurface = dd_cellline[ct2]['volumepsurface_final'] + ((dd_cellline['WT']['volumepsurface_final'] - dd_cellline[ct2]['volumepsurface_final']) * (time - ct2_kd_time)**hill_n) / (dd_cellline[ct2]['km_half_time']**hill_n + (time - ct1_kd_time)**hill_n)
                    # calculate surface
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
            # surface
            #cell.dict['volumepsurface'] = set_volumepsurface
            cell.targetSurface = int(set_surface)
            # adhesion
            #cell.dict['adhesion'] = set_adhesion
            self.adhesionFlexPlugin.setAdhesionMoleculeDensity(cell, 'cdh', set_adhesion)
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
            self.plot_win_adhesion.add_data_point("medium_adhesion", x=mcs, y=adhesion_init)
            self.plot_win_adhesion.add_data_point("wt_adhesion", x=mcs, y=adhesion_max)
            lambda_adhesion = self.get_xml_element("lamdaAdhesion").cdata
            self.plot_win_adhesion.add_data_point("lambda_adhesion", x=mcs, y=lambda_adhesion)
            # time
            self.plot_win_adhesion.add_data_point("experiment_time", x=mcs, y=time)
        # contact
        if (contact_verbose):
            # out
            self.plot_win_contact.add_data_point("total", x=mcs, y=sum(di_contact.values()))
            self.plot_win_contact.add_data_point("medium_ct1", x=mcs, y=di_contact['medium_ct1'] )
            self.plot_win_contact.add_data_point("medium_ct2", x=mcs, y=di_contact['medium_ct2'])
            self.plot_win_contact.add_data_point("ct1_ct1", x=mcs, y=di_contact['ct1_ct1'])
            self.plot_win_contact.add_data_point("ct1_ct2", x=mcs, y=di_contact['ct1_ct2'])
            self.plot_win_contact.add_data_point("ct2_ct2", x=mcs, y=int(di_contact['ct2_ct2']))
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
        Called after the last MCS to wrap up the simulation. Good place to close files and do post-processing
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
        pass

    def step(self, mcs):
        if (mcs==0):
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
                    set_generation_time = dd_cellline[ct1]['generation_time_final']
                else:
                    set_generation_time = dd_cellline[ct1]['generation_time_final'] + ((dd_cellline['WT']['generation_time_final'] - dd_cellline[ct1]['generation_time_final']) * (time - ct1_kd_time)**hill_n) / (dd_cellline[ct1]['km_half_time']**hill_n + (time - ct1_kd_time)**hill_n)

            #######
            # ct2 #
            #######
            elif (cell.type == self.CT2):
                # ct2_generation_time
                if (time < 0):
                    set_generation_time = dd_cellline[ct2]['generation_time_final']
                else:
                    set_generation_time = dd_cellline[ct2]['generation_time_final'] + ((dd_cellline['WT']['generation_time_final'] - dd_cellline[ct2]['generation_time_final']) * (time - ct2_kd_time)**hill_n) / (dd_cellline[ct2]['km_half_time']**hill_n + (time - ct2_kd_time)**hill_n)

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
            cell.dict['dc'] = 0  # set divison time counter to zero
            self.divide_cell_random_orientation(cell)

    def update_attributes(self):
        self.parent_cell.targetVolume /= 2.0  # reducing parent target volume
        self.clone_parent_2_child()

