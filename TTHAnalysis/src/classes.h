#include "CMGTools/TTHAnalysis/interface/SignedImpactParameter.h"
#include "CMGTools/TTHAnalysis/interface/DistributionRemapper.h"
#include "CMGTools/TTHAnalysis/interface/PdfWeightProducerTool.h"
#include "CMGTools/TTHAnalysis/interface/IgProfHook.h"
#include "CMGTools/TTHAnalysis/interface/CollectionSkimmer.h"
#include "CMGTools/TTHAnalysis/interface/CombinedObjectTags.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
namespace {
    struct dictionary {
        SignedImpactParameter sipc;
        DistributionRemapper remapper;
        PdfWeightProducerTool pdfw;
        SetupIgProfDumpHook hook;
        tensorflow::SessionOptions SessionOptions_;
    };
}
