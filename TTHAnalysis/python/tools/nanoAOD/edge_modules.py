MODULES = []
from CMGTools.TTHAnalysis.tools.nanoAOD.Edge_triggers import Triggers
MODULES.append( ('Triggers' , lambda : Triggers('Trigger','Filters')))

