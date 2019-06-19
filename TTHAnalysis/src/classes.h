#include "CMGTools/TTHAnalysis/interface/SignedImpactParameter.h"
#include "CMGTools/TTHAnalysis/interface/DistributionRemapper.h"
#include "CMGTools/TTHAnalysis/interface/PdfWeightProducerTool.h"
#include "CMGTools/TTHAnalysis/interface/IgProfHook.h"
#include "CMGTools/TTHAnalysis/interface/CollectionSkimmer.h"
#include "CMGTools/TTHAnalysis/interface/CombinedObjectTags.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#pragma GCC diagnostic pop
namespace {
    struct dictionary {
        SignedImpactParameter sipc;
        DistributionRemapper remapper;
        PdfWeightProducerTool pdfw;
        SetupIgProfDumpHook hook;
        tensorflow::SessionOptions SessionOptions_;
    };
}
