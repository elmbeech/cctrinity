
from cc3d import CompuCellSetup
        


from cctrinitySteppables import ConstraintInitializerSteppable

CompuCellSetup.register_steppable(steppable=ConstraintInitializerSteppable(frequency=1))




from cctrinitySteppables import GrowthSteppable

CompuCellSetup.register_steppable(steppable=GrowthSteppable(frequency=1))


CompuCellSetup.run()
