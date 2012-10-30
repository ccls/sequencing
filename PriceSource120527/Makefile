#specialParams = -fopenmp -g
specialParams = -fopenmp -O
#specialParams = -O
#specialParams = -fopenmp 

assemblerClasses = \
burrowsWheelerUtilities.o \
AssemblyException.o \
AlignmentScoreMatrix.o \
ParamsMinOverlap.o \
ParamsMinFractId.o \
ParamsAlignment.o \
ParamsMapping.o \
ParamsDeBruijn.o \
ParamsCycleManagement.o \
ScoredSeq.o \
ScoredSeqNested.o \
ScoredSeqPaired.o \
ScoredSeqMonoScore.o \
Read.o \
Alignment.o \
AlignmentNull.o \
AlignmentFull.o \
AlignmentUngapped.o \
AlignmentGapped.o \
DynamicProgrammingAligner.o \
ScoredSeqShallow.o \
ScoredSeqSubseq.o \
ScoredSeqNormalized.o \
ScoredSeqWithCarriers.o \
ScoredSeqFlip.o \
ReadFile.o \
ReadFileCommunicator.o \
ReadFileFastaSingle.o \
ReadFileFastqSingle.o \
ReadFilePriceqSingle.o \
ReadFileDual.o \
ReadFilePaired.o \
ReadFileIndex.o \
ParameterizedReadFile.o \
ParameterizedInitialFile.o \
EcoFilter.o \
RepeatDetector.o \
ReadPairFilter.o \
OutputFile.o \
OutputFileNull.o \
OutputFileFasta.o \
OutputFilePriceq.o \
OutputFileMulti.o \
MatchSeqTest.o \
MatchOffsetTest.o \
MatchAlignTest.o \
ScoredSeqCollectionBwt.o \
AssemblyJob.o \
AssemblyJobFactory.o \
AssemblyJobIterativeBwt.o \
AssemblyJobNull.o \
AssemblyJobNullCopy.o \
AssemblyJobHierarchy.o \
AssemblyJobPerfect.o \
AssemblyJobNearPerfect.o \
AssemblyJobSubset.o \
AssemblyJobSsDeBruijn.o \
AssemblyGraphNode.o \
AssemblyJobGraph.o \
MapReadSetupHelper.o \
ExtendJobMapper.o \
ExtendJobCreator.o \
ExtendCycle.o \
Assembler.o \
AssemblerListener.o \
AssemblerListenerNull.o \
AssemblerListenerStd.o \
AssemblerListenerVerbose.o \
AssemblerListenerCarrier.o \
# \
nucleicAcidUtilities.o \
AssemblyJobPair.o \
AssemblerListenerFile.o \
ScoredSeqCollection.o \
AssemblyJobIterative.o \
AssemblyJobRedundancy.o \
ScoredSeqAlignCarriers.o \
AssemblyJobHelperStrategies.o \


PriceTI: $(assemblerClasses) PriceTI.o
	g++ $(specialParams) -o PriceTI PriceTI.o $(assemblerClasses)

PriceTI.o: PriceTI.cpp
	g++ $(specialParams) -c PriceTI.cpp

# GENERAL USE
burrowsWheelerUtilities.o: burrowsWheelerUtilities.h burrowsWheelerUtilities.cpp
	g++ $(specialParams) -c burrowsWheelerUtilities.cpp -o burrowsWheelerUtilities.o

AssemblyException.o: AssemblyException.h AssemblyException.cpp
	g++ $(specialParams) -c AssemblyException.cpp -o AssemblyException.o

AlignmentScoreMatrix.o: AlignmentScoreMatrix.h AlignmentScoreMatrix.cpp
	g++ $(specialParams) -c AlignmentScoreMatrix.cpp -o AlignmentScoreMatrix.o

ParamsMinOverlap.o: ParamsMinOverlap.h ParamsMinOverlap.cpp
	g++ $(specialParams) -c ParamsMinOverlap.cpp -o ParamsMinOverlap.o

ParamsMinFractId.o: ParamsMinFractId.h ParamsMinFractId.cpp
	g++ $(specialParams) -c ParamsMinFractId.cpp -o ParamsMinFractId.o

ParamsAlignment.o: ParamsAlignment.h ParamsAlignment.cpp
	g++ $(specialParams) -c ParamsAlignment.cpp -o ParamsAlignment.o

ParamsMapping.o: ParamsMapping.h ParamsMapping.cpp
	g++ $(specialParams) -c ParamsMapping.cpp -o ParamsMapping.o

ParamsDeBruijn.o: ParamsDeBruijn.h ParamsDeBruijn.cpp
	g++ $(specialParams) -c ParamsDeBruijn.cpp -o ParamsDeBruijn.o

EcoFilter.o: EcoFilter.h EcoFilter.cpp
	g++ $(specialParams) -c EcoFilter.cpp -o EcoFilter.o

ReadPairFilter.o: ReadPairFilter.h ReadPairFilter.cpp
	g++ $(specialParams) -c ReadPairFilter.cpp -o ReadPairFilter.o

ParamsCycleManagement.o: ParamsCycleManagement.h ParamsCycleManagement.cpp
	g++ $(specialParams) -c ParamsCycleManagement.cpp -o ParamsCycleManagement.o

MatchSeqTest.o: MatchSeqTest.h MatchSeqTest.cpp
	g++ $(specialParams) -c MatchSeqTest.cpp -o MatchSeqTest.o

MatchOffsetTest.o: MatchOffsetTest.h MatchOffsetTest.cpp
	g++ $(specialParams) -c MatchOffsetTest.cpp -o MatchOffsetTest.o

MatchAlignTest.o: MatchAlignTest.h MatchAlignTest.cpp
	g++ $(specialParams) -c MatchAlignTest.cpp -o MatchAlignTest.o

RepeatDetector.o: RepeatDetector.h RepeatDetector.cpp
	g++ $(specialParams) -c RepeatDetector.cpp -o RepeatDetector.o



# HIGER-ORDER OPERATIONS
ScoredSeqCollection.o: ScoredSeqCollection.h ScoredSeqCollection.cpp
	g++ $(specialParams) -c ScoredSeqCollection.cpp -o ScoredSeqCollection.o

ScoredSeqCollectionBwt.o: ScoredSeqCollectionBwt.h ScoredSeqCollectionBwt.cpp
	g++ $(specialParams) -c ScoredSeqCollectionBwt.cpp -o ScoredSeqCollectionBwt.o

AssemblyJob.o: AssemblyJob.h AssemblyJob.cpp
	g++ $(specialParams) -c AssemblyJob.cpp -o AssemblyJob.o

AssemblyJobFactory.o: AssemblyJobFactory.h AssemblyJobFactory.cpp
	g++ $(specialParams) -c AssemblyJobFactory.cpp -o AssemblyJobFactory.o

AssemblyJobIterative.o: AssemblyJobIterative.h AssemblyJobIterative.cpp
	g++ $(specialParams) -c AssemblyJobIterative.cpp -o AssemblyJobIterative.o

AssemblyJobIterativeBwt.o: AssemblyJobIterativeBwt.h AssemblyJobIterativeBwt.cpp
	g++ $(specialParams) -c AssemblyJobIterativeBwt.cpp -o AssemblyJobIterativeBwt.o

AssemblyJobNull.o: AssemblyJobNull.h AssemblyJobNull.cpp
	g++ $(specialParams) -c AssemblyJobNull.cpp -o AssemblyJobNull.o

AssemblyJobNullCopy.o: AssemblyJobNullCopy.h AssemblyJobNullCopy.cpp
	g++ $(specialParams) -c AssemblyJobNullCopy.cpp -o AssemblyJobNullCopy.o

AssemblyGraphNode.o: AssemblyGraphNode.h AssemblyGraphNode.cpp
	g++ $(specialParams) -c AssemblyGraphNode.cpp -o AssemblyGraphNode.o

AssemblyJobGraph.o: AssemblyJobGraph.h AssemblyJobGraph.cpp
	g++ $(specialParams) -c AssemblyJobGraph.cpp -o AssemblyJobGraph.o

AssemblyJobHierarchy.o: AssemblyJobHierarchy.h AssemblyJobHierarchy.cpp
	g++ $(specialParams) -c AssemblyJobHierarchy.cpp -o AssemblyJobHierarchy.o

AssemblyJobRedundancy.o: AssemblyJobRedundancy.h AssemblyJobRedundancy.cpp
	g++ $(specialParams) -c AssemblyJobRedundancy.cpp -o AssemblyJobRedundancy.o

AssemblyJobPerfect.o: AssemblyJobPerfect.h AssemblyJobPerfect.cpp
	g++ $(specialParams) -c AssemblyJobPerfect.cpp -o AssemblyJobPerfect.o

AssemblyJobNearPerfect.o: AssemblyJobNearPerfect.h AssemblyJobNearPerfect.cpp
	g++ $(specialParams) -c AssemblyJobNearPerfect.cpp -o AssemblyJobNearPerfect.o

AssemblyJobSubset.o: AssemblyJobSubset.h AssemblyJobSubset.cpp
	g++ $(specialParams) -c AssemblyJobSubset.cpp -o AssemblyJobSubset.o

AssemblyJobSsDeBruijn.o: AssemblyJobSsDeBruijn.h AssemblyJobSsDeBruijn.cpp
	g++ $(specialParams) -c AssemblyJobSsDeBruijn.cpp -o AssemblyJobSsDeBruijn.o

AssemblyJobHelperStrategies.o: AssemblyJobHelperStrategies.h AssemblyJobHelperStrategies.cpp
	g++ $(specialParams) -c AssemblyJobHelperStrategies.cpp -o AssemblyJobHelperStrategies.o


# SCORED SEQS
ScoredSeq.o: ScoredSeq.h ScoredSeq.cpp
	g++ $(specialParams) -c ScoredSeq.cpp -o ScoredSeq.o

ScoredSeqNested.o: ScoredSeqNested.h ScoredSeqNested.cpp
	g++ $(specialParams) -c ScoredSeqNested.cpp -o ScoredSeqNested.o

ScoredSeqPaired.o: ScoredSeqPaired.h ScoredSeqPaired.cpp
	g++ $(specialParams) -c ScoredSeqPaired.cpp -o ScoredSeqPaired.o

ScoredSeqShallow.o: ScoredSeqShallow.h ScoredSeqShallow.cpp
	g++ $(specialParams) -c ScoredSeqShallow.cpp -o ScoredSeqShallow.o

ScoredSeqSubseq.o: ScoredSeqSubseq.h ScoredSeqSubseq.cpp
	g++ $(specialParams) -c ScoredSeqSubseq.cpp -o ScoredSeqSubseq.o

ScoredSeqMonoScore.o: ScoredSeqMonoScore.h ScoredSeqMonoScore.cpp
	g++ $(specialParams) -c ScoredSeqMonoScore.cpp -o ScoredSeqMonoScore.o

ScoredSeqNormalized.o: ScoredSeqNormalized.h ScoredSeqNormalized.cpp
	g++ $(specialParams) -c ScoredSeqNormalized.cpp -o ScoredSeqNormalized.o

ScoredSeqWithCarriers.o: ScoredSeqWithCarriers.h ScoredSeqWithCarriers.cpp
	g++ $(specialParams) -c ScoredSeqWithCarriers.cpp -o ScoredSeqWithCarriers.o

ScoredSeqAlignCarriers.o: ScoredSeqAlignCarriers.h ScoredSeqAlignCarriers.cpp
	g++ $(specialParams) -c ScoredSeqAlignCarriers.cpp -o ScoredSeqAlignCarriers.o

ScoredSeqFlip.o: ScoredSeqFlip.h ScoredSeqFlip.cpp
	g++ $(specialParams) -c ScoredSeqFlip.cpp -o ScoredSeqFlip.o

Read.o: Read.h Read.cpp
	g++ $(specialParams) -c Read.cpp -o Read.o


# ALIGNMENTS
Alignment.o: Alignment.h Alignment.cpp
	g++ $(specialParams) -c Alignment.cpp -o Alignment.o

AlignmentNull.o: AlignmentNull.h AlignmentNull.cpp
	g++ $(specialParams) -c AlignmentNull.cpp -o AlignmentNull.o

AlignmentFull.o: AlignmentFull.h AlignmentFull.cpp
	g++ $(specialParams) -c AlignmentFull.cpp -o AlignmentFull.o

AlignmentUngapped.o: AlignmentUngapped.h AlignmentUngapped.cpp
	g++ $(specialParams) -c AlignmentUngapped.cpp -o AlignmentUngapped.o

AlignmentGapped.o: AlignmentGapped.h AlignmentGapped.cpp
	g++ $(specialParams) -c AlignmentGapped.cpp -o AlignmentGapped.o

DynamicProgrammingAligner.o: DynamicProgrammingAligner.h DynamicProgrammingAligner.cpp
	g++ $(specialParams) -c DynamicProgrammingAligner.cpp -o DynamicProgrammingAligner.o


# READ FILES
ReadFile.o: ReadFile.h ReadFile.cpp
	g++ $(specialParams) -c ReadFile.cpp -o ReadFile.o

ReadFileCommunicator.o: ReadFileCommunicator.h ReadFileCommunicator.cpp
	g++ $(specialParams) -c ReadFileCommunicator.cpp -o ReadFileCommunicator.o

ReadFileFastqSingle.o: ReadFileFastqSingle.h ReadFileFastqSingle.cpp
	g++ $(specialParams) -c ReadFileFastqSingle.cpp -o ReadFileFastqSingle.o

ReadFilePriceqSingle.o: ReadFilePriceqSingle.h ReadFilePriceqSingle.cpp
	g++ $(specialParams) -c ReadFilePriceqSingle.cpp -o ReadFilePriceqSingle.o

ReadFileFastaSingle.o: ReadFileFastaSingle.h ReadFileFastaSingle.cpp
	g++ $(specialParams) -c ReadFileFastaSingle.cpp -o ReadFileFastaSingle.o

ReadFileDual.o: ReadFileDual.h ReadFileDual.cpp
	g++ $(specialParams) -c ReadFileDual.cpp -o ReadFileDual.o

ReadFilePaired.o: ReadFilePaired.h ReadFilePaired.cpp
	g++ $(specialParams) -c ReadFilePaired.cpp -o ReadFilePaired.o

ReadFileIndex.o: ReadFileIndex.h ReadFileIndex.cpp
	g++ $(specialParams) -c ReadFileIndex.cpp -o ReadFileIndex.o

OutputFile.o: OutputFile.h OutputFile.cpp
	g++ $(specialParams) -c OutputFile.cpp -o OutputFile.o

OutputFileNull.o: OutputFileNull.h OutputFileNull.cpp
	g++ $(specialParams) -c OutputFileNull.cpp -o OutputFileNull.o

OutputFileFasta.o: OutputFileFasta.h OutputFileFasta.cpp
	g++ $(specialParams) -c OutputFileFasta.cpp -o OutputFileFasta.o

OutputFilePriceq.o: OutputFilePriceq.h OutputFilePriceq.cpp
	g++ $(specialParams) -c OutputFilePriceq.cpp -o OutputFilePriceq.o

OutputFileMulti.o: OutputFileMulti.h OutputFileMulti.cpp
	g++ $(specialParams) -c OutputFileMulti.cpp -o OutputFileMulti.o

ParameterizedInitialFile.o: ParameterizedInitialFile.h ParameterizedInitialFile.cpp
	g++ $(specialParams) -c ParameterizedInitialFile.cpp -o ParameterizedInitialFile.o

ParameterizedReadFile.o: ParameterizedReadFile.h ParameterizedReadFile.cpp
	g++ $(specialParams) -c ParameterizedReadFile.cpp -o ParameterizedReadFile.o

#ReadFileFactory.o: ReadFileFactory.h ReadFileFactory.cpp
#	g++ -c ReadFileFactory.cpp -o ReadFileFactory.o


# HIGHER-ORDER OPERATIONS
MapReadSetupHelper.o: MapReadSetupHelper.h MapReadSetupHelper.cpp
	g++ $(specialParams) -c MapReadSetupHelper.cpp -o MapReadSetupHelper.o

ExtendJobMapper.o: ExtendJobMapper.h ExtendJobMapper.cpp
	g++ $(specialParams) -c ExtendJobMapper.cpp -o ExtendJobMapper.o

ExtendJobCreator.o: ExtendJobCreator.h ExtendJobCreator.cpp
	g++ $(specialParams) -c ExtendJobCreator.cpp -o ExtendJobCreator.o

ExtendCycle.o: ExtendCycle.h ExtendCycle.cpp
	g++ $(specialParams) -c ExtendCycle.cpp -o ExtendCycle.o

Assembler.o: Assembler.h Assembler.cpp
	g++ $(specialParams) -c Assembler.cpp -o Assembler.o



# COMMUNICATIONS
AssemblerListener.o: AssemblerListener.h AssemblerListener.cpp
	g++ $(specialParams) -c AssemblerListener.cpp -o AssemblerListener.o

AssemblerListenerNull.o: AssemblerListenerNull.h AssemblerListenerNull.cpp
	g++ $(specialParams) -c AssemblerListenerNull.cpp -o AssemblerListenerNull.o

AssemblerListenerStd.o: AssemblerListenerStd.h AssemblerListenerStd.cpp
	g++ $(specialParams) -c AssemblerListenerStd.cpp -o AssemblerListenerStd.o

AssemblerListenerVerbose.o: AssemblerListenerVerbose.h AssemblerListenerVerbose.cpp
	g++ $(specialParams) -c AssemblerListenerVerbose.cpp -o AssemblerListenerVerbose.o

AssemblerListenerCarrier.o: AssemblerListenerCarrier.h AssemblerListenerCarrier.cpp
	g++ $(specialParams) -c AssemblerListenerCarrier.cpp -o AssemblerListenerCarrier.o


clean:
	rm -f $(assemblerClasses) PriceTI.o TrialMain3.o

