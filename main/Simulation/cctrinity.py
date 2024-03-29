# library
from cc3d import CompuCellSetup
from cctrinitySteppables import GrowthSteppable, MitosisSteppable

# register steppable
CompuCellSetup.register_steppable(steppable=GrowthSteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))

# run cc3d
CompuCellSetup.run()
