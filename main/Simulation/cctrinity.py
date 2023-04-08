# library
from cc3d import CompuCellSetup
from cctrinitySteppables import ConstraintInitializerSteppable
from cctrinitySteppables import GrowthSteppable
from cctrinitySteppables import MitosisSteppable

# register steppable
CompuCellSetup.register_steppable(steppable=ConstraintInitializerSteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=GrowthSteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))

# run cc3d
CompuCellSetup.run()
