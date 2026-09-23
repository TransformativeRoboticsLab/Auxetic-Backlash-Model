import Std

namespace RadDiscontinuities

/-!
Variational and contact-mechanics scaffold.

This file intentionally separates proved toy energy facts from ambitious
real-valued mechanics targets. The `Nat` energies provide small compiled
nonnegativity theorems. The generic scalar structures record the premises a
future real-valued/mathlib-backed model must satisfy before we treat RAD contact,
pin-hole clearance, hinges, or springs as proved physics.
-/

def springEnergyNat (stiffness displacement : Nat) : Nat :=
  stiffness * displacement * displacement

def hingeEnergyNat (stiffness angleError : Nat) : Nat :=
  stiffness * angleError * angleError

def contactPenaltyNat (stiffness penetration : Nat) : Nat :=
  stiffness * penetration * penetration

def clearanceExcessNat (clearance command : Nat) : Nat :=
  command - clearance

def verticalResidualStepNat (clearance gain command : Nat) : Nat :=
  gain * clearanceExcessNat clearance command

def contactPenaltyFromClearanceNat
    (stiffness clearance command : Nat) : Nat :=
  contactPenaltyNat stiffness (clearanceExcessNat clearance command)

def loadWorkMagnitudeNat (load displacement : Nat) : Nat :=
  load * displacement

def mechanicalStoredEnergyNat
    (spring hinge lock contact loadMagnitude : Nat) : Nat :=
  spring + hinge + lock + contact + loadMagnitude

def absoluteErrorNat (measured simulated : Nat) : Nat :=
  if simulated <= measured then measured - simulated else simulated - measured

def measuredWorkResidualNat (measured simulated : Nat) : Nat :=
  absoluteErrorNat measured simulated

def signedLoadWorkInt (load displacement : Int) : Int :=
  -load * displacement

def absoluteErrorInt (measured simulated : Int) : Nat :=
  (measured - simulated).natAbs

def signedWorkResidualInt (measured simulated : Int) : Nat :=
  absoluteErrorInt measured simulated

def measuredEnergyResidualTripleNat
    (measuredSigned simulatedSigned
      measuredMagnitude simulatedMagnitude
      measuredContact simulatedContact : Nat) : Nat :=
  measuredWorkResidualNat measuredSigned simulatedSigned
    + measuredWorkResidualNat measuredMagnitude simulatedMagnitude
    + measuredWorkResidualNat measuredContact simulatedContact

def verticalLoadScenarioPassNat
    (missingMeasurements signedError magnitudeError contactError tolerance : Nat) : Bool :=
  missingMeasurements == 0
    && signedError <= tolerance
    && magnitudeError <= tolerance
    && contactError <= tolerance

def signedEnergyResidualTripleInt
    (measuredSigned simulatedSigned
      measuredMagnitude simulatedMagnitude
      measuredContact simulatedContact : Int) : Nat :=
  signedWorkResidualInt measuredSigned simulatedSigned
    + signedWorkResidualInt measuredMagnitude simulatedMagnitude
    + signedWorkResidualInt measuredContact simulatedContact

structure ScaledNatQuantity where
  numerator : Nat
  denominator : Nat

structure ScaledIntQuantity where
  numerator : Int
  denominator : Nat

structure MeasurementUnitScale where
  numerator : Nat
  denominator : Nat

structure HardwareProfileCoverage where
  measuredFields : Nat
  totalFields : Nat

def hardwareProfileCoverageComplete (profile : HardwareProfileCoverage) : Bool :=
  profile.measuredFields == profile.totalFields

def hardwareProfileCoverageMissing (profile : HardwareProfileCoverage) : Nat :=
  profile.totalFields - profile.measuredFields

structure CalibrationFitResidualNat where
  sampleCount : Nat
  rawResidual : Nat
  fittedResidual : Nat
  tolerance : Nat

def calibrationFitResidualPassNat (fit : CalibrationFitResidualNat) : Prop :=
  fit.sampleCount > 0 ∧ fit.fittedResidual <= fit.tolerance

structure CalibrationModelProfileUpdateNat where
  sampleCount : Nat
  lowerBound : Nat
  upperBound : Nat
  proposed : Nat

structure CalibrationModelProfileCandidateNat where
  beforeResidual : Nat
  afterResidual : Nat
  missingBefore : Nat
  missingAfter : Nat
  appliedUpdates : Nat

structure CalibrationModelProfileHoldoutNat where
  fitBeforeResidual : Nat
  fitAfterResidual : Nat
  holdoutBeforeResidual : Nat
  holdoutAfterResidual : Nat
  fitMissingBefore : Nat
  fitMissingAfter : Nat
  holdoutMissingBefore : Nat
  holdoutMissingAfter : Nat
  appliedUpdates : Nat

structure CalibrationTrainHoldoutSplitNat where
  fitSamples : Nat
  holdoutSamples : Nat
  overlapSamples : Nat
  frozenProfile : Nat

structure CalibrationTrainHoldoutProvenanceNat where
  fitDatasetIdKnown : Nat
  holdoutDatasetIdKnown : Nat
  distinctDatasetIds : Nat
  distinctSourceFileIds : Nat
  fitRoleMarked : Nat
  holdoutRoleMarked : Nat
  profileFrozen : Nat
  holdoutProfileMatchesFrozenProfile : Nat

structure CalibrationBenchProtocolCoverageNat where
  scenarioCount : Nat
  instrumentCount : Nat
  outputArtifactCount : Nat
  measurementColumnCount : Nat
  datasetRoleCount : Nat

structure CalibrationBenchPacketCompletenessNat where
  notebookArtifacts : Nat
  protocolArtifacts : Nat
  fitTemplateArtifacts : Nat
  holdoutTemplateArtifacts : Nat
  validationInstructionArtifacts : Nat
  manifestEntries : Nat

structure CalibrationBenchExecutedValidationNat where
  fitSamples : Nat
  holdoutSamples : Nat
  appliedUpdates : Nat
  residualValidationPass : Nat
  independentValidationPass : Nat
  missingEvidence : Nat

structure PhysicalValidationReadinessNat where
  calibrationExecutionReady : Nat
  verticalLoadScenarios : Nat
  verticalLoadComparisonPass : Nat
  loadProxyTerms : Nat
  contactProxyTerms : Nat
  clearanceConfigured : Nat
  missingEvidence : Nat

structure ContactStateAbstractionNat where
  bodyCount : Nat
  pinCount : Nat
  holeCount : Nat
  clearancePairCount : Nat
  contactStateRecords : Nat
  penaltyTerms : Nat
  clearanceConfigured : Nat
  missingEvidence : Nat

structure ContactGraphConsistencyNat where
  activeBodies : Nat
  contactRecords : Nat
  activeEdges : Nat
  removedIncidentActiveEdges : Nat
  removedActiveContacts : Nat
  groupSupportCells : Nat
  groupSupportRecords : Nat
  missingEvidence : Nat

structure PhysicalRealizationMapNat where
  abstractOperators : Nat
  realizedOperators : Nat
  supportRecords : Nat
  stateEffectRecords : Nat
  claimLabels : Nat
  contactGraphReady : Nat
  missingEvidence : Nat

structure ExternalPhysicsEngineAuditNat where
  engineCandidates : Nat
  availableCandidates : Nat
  requiredFeatureRecords : Nat
  scenarioRecords : Nat
  contactModelRecords : Nat
  independentToolRecords : Nat
  missingEvidence : Nat

structure ExternalPhysicsModelExportNat where
  modelRecords : Nat
  bodyRecords : Nat
  fixedBodyRecords : Nat
  gravityRecords : Nat
  xmlBytes : Nat
  missingEvidence : Nat

structure ExternalContactGeometryNat where
  pinRecords : Nat
  holeRecords : Nat
  clearanceRecords : Nat
  contactPairRecords : Nat
  activeContactPairs : Nat
  xmlBytes : Nat
  missingEvidence : Nat

structure ExternalContactParameterProfileNat where
  contactPairRecords : Nat
  parameterRecords : Nat
  frictionRecords : Nat
  solverParameterRecords : Nat
  stiffnessRecords : Nat
  dampingRecords : Nat
  xmlAttributeBytes : Nat
  missingEvidence : Nat

structure ContactParameterCalibrationPacketNat where
  contactPairRecords : Nat
  parameterRecords : Nat
  measurementColumns : Nat
  fitTemplateRows : Nat
  holdoutTemplateRows : Nat
  artifactManifestEntries : Nat
  profileEvidenceReady : Nat
  missingEvidence : Nat

structure ContactParameterBenchValidationNat where
  packetReady : Nat
  fitRows : Nat
  holdoutRows : Nat
  completedMeasurements : Nat
  parameterResiduals : Nat
  fitPass : Nat
  holdoutPass : Nat
  independentHoldoutPass : Nat
  missingEvidence : Nat

structure ContactParameterIntervalCalibrationNat where
  benchValidationReady : Nat
  parameterIntervals : Nat
  acceptedParameterIntervals : Nat
  uncertaintyRecords : Nat
  holdoutAgreementRecords : Nat
  simulatorParametersInsideBounds : Nat
  missingEvidence : Nat

structure ExternalPhysicsRunNat where
  engineAvailable : Nat
  modelExportReady : Nat
  bodyRecords : Nat
  resultBodyRecords : Nat
  stepsRequested : Nat
  stepsCompleted : Nat
  missingEvidence : Nat

structure ExternalPhysicsComparisonNat where
  runReady : Nat
  simulatorRecords : Nat
  externalRecords : Nat
  matchedRecords : Nat
  toleranceRecords : Nat
  withinTolerance : Nat
  missingEvidence : Nat

structure EquilibriumRelationNat where
  realizationReady : Nat
  solverSuccess : Nat
  contactStateReady : Nat
  energyTerms : Nat
  nonnegativeEnergyTerms : Nat
  balanceResidual : Nat
  residualTolerance : Nat
  missingEvidence : Nat

structure ReachableEquilibriumControllabilityNat where
  equilibriumReady : Nat
  actuatorBasis : Nat
  responseColumns : Nat
  targetCells : Nat
  reachableResponses : Nat
  topologyEvidence : Nat
  fullTargetRequired : Nat
  targetUnderactuated : Nat
  missingEvidence : Nat

structure ReachableEquilibriumBenchProtocolNat where
  reachabilityReady : Nat
  scenarioCount : Nat
  actuatorColumnTrials : Nat
  targetObservationCells : Nat
  measurementColumns : Nat
  topologyControlPolicy : Nat
  passFailCriteria : Nat
  missingEvidence : Nat

structure ReachableEquilibriumBenchValidationNat where
  protocolReady : Nat
  resultRows : Nat
  targetMeasurements : Nat
  completedMeasurements : Nat
  reachabilityChecks : Nat
  topologyLeakageChecks : Nat
  groupSequenceChecks : Nat
  comparisonPass : Nat
  missingEvidence : Nat

structure ReachableEquilibriumAmplitudeCalibrationNat where
  benchValidationReady : Nat
  repeatedTrialGroups : Nat
  amplitudeEstimates : Nat
  residualFieldCells : Nat
  uncertaintyBands : Nat
  topologyLeakageBands : Nat
  groupSequenceResiduals : Nat
  missingEvidence : Nat

structure ReachableEquilibriumEmpiricalProfileNat where
  amplitudeCalibrationReady : Nat
  boundedProposals : Nat
  safeProposals : Nat
  toleranceRecords : Nat
  uncertaintyRecords : Nat
  holdoutHooks : Nat
  missingEvidence : Nat

structure ReachableEquilibriumProfileInverseNat where
  empiricalProfileReady : Nat
  safeProfileProposals : Nat
  targetCells : Nat
  inverseSolveSuccess : Nat
  residualRecords : Nat
  scoreTerms : Nat
  readOnlyProfileUse : Nat
  missingEvidence : Nat

structure ReachableEquilibriumProfileInverseAcceptanceNat where
  profileInverseReady : Nat
  scoreRecords : Nat
  residualScore : Nat
  residualLimit : Nat
  bandFailures : Nat
  maxBandFailures : Nat
  actuatorCount : Nat
  actuatorLimit : Nat
  acceptedForPreview : Nat
  missingEvidence : Nat

structure ReachableEquilibriumProfileInversePreviewPacketNat where
  acceptanceReady : Nat
  commandRecords : Nat
  eventRecords : Nat
  targetRecords : Nat
  residualRecords : Nat
  readOnlyPacket : Nat
  missingEvidence : Nat

structure ReachableEquilibriumProfileInversePreviewReplayNat where
  packetReady : Nat
  commandRecords : Nat
  replayedCommands : Nat
  eventCoverage : Nat
  simulationRecords : Nat
  targetResidualRecords : Nat
  residualAgreement : Nat
  readOnlyReplay : Nat
  missingEvidence : Nat

structure ReachableEquilibriumProfileInversePreviewPhysicalNat where
  replayReady : Nat
  physicalSolverSuccess : Nat
  commandRecords : Nat
  targetResidualRecords : Nat
  modelComparisonRecords : Nat
  energyRecords : Nat
  readOnlyPhysicalPreview : Nat
  missingEvidence : Nat

def calibrationModelProfileUpdateSafeNat
    (update : CalibrationModelProfileUpdateNat) : Prop :=
  update.sampleCount > 0 ∧
    update.lowerBound <= update.proposed ∧
    update.proposed <= update.upperBound

def calibrationModelProfileCandidateImprovesNat
    (candidate : CalibrationModelProfileCandidateNat) : Prop :=
  And
    (candidate.appliedUpdates > 0)
    (And
      (candidate.afterResidual <= candidate.beforeResidual)
      (candidate.missingAfter <= candidate.missingBefore))

def calibrationModelProfileCandidateScoreNat
    (candidate : CalibrationModelProfileCandidateNat) : Nat :=
  candidate.afterResidual

def calibrationModelProfileHoldoutPassNat
    (validation : CalibrationModelProfileHoldoutNat) : Prop :=
  And
    (validation.appliedUpdates > 0)
    (And
      (validation.fitAfterResidual <= validation.fitBeforeResidual)
      (And
        (validation.holdoutAfterResidual <= validation.holdoutBeforeResidual)
        (And
          (validation.fitMissingAfter <= validation.fitMissingBefore)
          (validation.holdoutMissingAfter <= validation.holdoutMissingBefore))))

def calibrationModelProfileHoldoutScoreNat
    (validation : CalibrationModelProfileHoldoutNat) : Nat :=
  validation.holdoutAfterResidual

def calibrationTrainHoldoutSplitReadyNat
    (split : CalibrationTrainHoldoutSplitNat) : Prop :=
  And
    (split.fitSamples > 0)
    (And
      (split.holdoutSamples > 0)
      (And
        (split.overlapSamples = 0)
        (split.frozenProfile > 0)))

def calibrationTrainHoldoutProvenanceReadyNat
    (provenance : CalibrationTrainHoldoutProvenanceNat) : Prop :=
  And
    (provenance.fitDatasetIdKnown > 0)
    (And
      (provenance.holdoutDatasetIdKnown > 0)
      (And
        (provenance.distinctDatasetIds > 0)
        (And
          (provenance.distinctSourceFileIds > 0)
          (And
            (provenance.fitRoleMarked > 0)
            (And
              (provenance.holdoutRoleMarked > 0)
              (And
                (provenance.profileFrozen > 0)
                (provenance.holdoutProfileMatchesFrozenProfile > 0)))))))

def calibrationBenchProtocolCoverageReadyNat
    (coverage : CalibrationBenchProtocolCoverageNat) : Prop :=
  And
    (coverage.scenarioCount > 0)
    (And
      (coverage.instrumentCount > 0)
      (And
        (coverage.outputArtifactCount > 0)
        (And
          (coverage.measurementColumnCount > 0)
          (coverage.datasetRoleCount >= 2))))

def calibrationBenchPacketCompleteNat
    (packet : CalibrationBenchPacketCompletenessNat) : Prop :=
  And
    (packet.notebookArtifacts > 0)
    (And
      (packet.protocolArtifacts > 0)
      (And
        (packet.fitTemplateArtifacts > 0)
        (And
          (packet.holdoutTemplateArtifacts > 0)
          (And
            (packet.validationInstructionArtifacts > 0)
            (packet.manifestEntries >= 5)))))

def calibrationBenchExecutedValidationReadyNat
    (validation : CalibrationBenchExecutedValidationNat) : Prop :=
  And
    (validation.fitSamples > 0)
    (And
      (validation.holdoutSamples > 0)
      (And
        (validation.appliedUpdates > 0)
        (And
          (validation.residualValidationPass > 0)
          (And
            (validation.independentValidationPass > 0)
            (validation.missingEvidence = 0)))))

def physicalValidationReadyNat
    (readiness : PhysicalValidationReadinessNat) : Prop :=
  And
    (readiness.calibrationExecutionReady > 0)
    (And
      (readiness.verticalLoadScenarios > 0)
      (And
        (readiness.verticalLoadComparisonPass > 0)
        (And
          (readiness.loadProxyTerms > 0)
          (And
            (readiness.contactProxyTerms > 0)
            (And
              (readiness.clearanceConfigured > 0)
              (readiness.missingEvidence = 0))))))

def contactStateAbstractionReadyNat
    (contact : ContactStateAbstractionNat) : Prop :=
  contact.bodyCount > 0 ∧
    contact.pinCount = contact.bodyCount ∧
    contact.holeCount = contact.bodyCount ∧
    contact.clearancePairCount = contact.bodyCount ∧
    contact.contactStateRecords >= contact.bodyCount ∧
    contact.penaltyTerms >= contact.bodyCount ∧
    contact.clearanceConfigured > 0 ∧
    contact.missingEvidence = 0

def contactGraphConsistentNat
    (graph : ContactGraphConsistencyNat) : Prop :=
  And
    (graph.activeBodies > 0)
    (And
      (graph.contactRecords >= graph.activeBodies)
      (And
        (graph.removedIncidentActiveEdges = 0)
        (And
          (graph.removedActiveContacts = 0)
          (And
            (graph.groupSupportRecords >= graph.groupSupportCells)
            (graph.missingEvidence = 0)))))

def physicalRealizationMapReadyNat
    (realization : PhysicalRealizationMapNat) : Prop :=
  And
    (realization.abstractOperators > 0)
    (And
      (realization.realizedOperators >= realization.abstractOperators)
      (And
        (realization.supportRecords >= realization.abstractOperators)
        (And
          (realization.stateEffectRecords >= realization.abstractOperators)
          (And
            (realization.claimLabels >= realization.abstractOperators)
            (And
              (realization.contactGraphReady > 0)
              (realization.missingEvidence = 0))))))

def externalPhysicsEngineAuditReadyNat
    (audit : ExternalPhysicsEngineAuditNat) : Prop :=
  And
    (audit.engineCandidates > 0)
    (And
      (audit.availableCandidates > 0)
      (And
        (audit.requiredFeatureRecords > 0)
        (And
          (audit.scenarioRecords > 0)
          (And
            (audit.contactModelRecords > 0)
            (And
              (audit.independentToolRecords > 0)
              (audit.missingEvidence = 0))))))

def externalPhysicsModelExportReadyNat
    (modelExport : ExternalPhysicsModelExportNat) : Prop :=
  And
    (modelExport.modelRecords > 0)
    (And
      (modelExport.bodyRecords > 0)
      (And
        (modelExport.fixedBodyRecords > 0)
        (And
          (modelExport.gravityRecords > 0)
          (And
            (modelExport.xmlBytes > 0)
            (modelExport.missingEvidence = 0)))))

def externalContactGeometryReadyNat
    (geometry : ExternalContactGeometryNat) : Prop :=
  And
    (geometry.pinRecords > 0)
    (And
      (geometry.holeRecords >= geometry.pinRecords)
      (And
        (geometry.clearanceRecords >= geometry.pinRecords)
        (And
          (geometry.contactPairRecords >= geometry.pinRecords)
          (And
            (geometry.activeContactPairs >= geometry.pinRecords)
            (And
              (geometry.xmlBytes > 0)
              (geometry.missingEvidence = 0))))))

def externalContactParameterProfileReadyNat
    (profile : ExternalContactParameterProfileNat) : Prop :=
  And
    (profile.contactPairRecords > 0)
    (And
      (profile.parameterRecords >= profile.contactPairRecords)
      (And
        (profile.frictionRecords >= profile.contactPairRecords)
        (And
          (profile.solverParameterRecords >= profile.contactPairRecords)
          (And
            (profile.stiffnessRecords >= profile.contactPairRecords)
            (And
              (profile.dampingRecords >= profile.contactPairRecords)
              (And
                (profile.xmlAttributeBytes > 0)
                (profile.missingEvidence = 0)))))))

def contactParameterCalibrationPacketCompleteNat
    (packet : ContactParameterCalibrationPacketNat) : Prop :=
  And
    (packet.contactPairRecords > 0)
    (And
      (packet.parameterRecords >= packet.contactPairRecords)
      (And
        (packet.measurementColumns > 0)
        (And
          (packet.fitTemplateRows >= packet.contactPairRecords)
          (And
            (packet.holdoutTemplateRows >= packet.contactPairRecords)
            (And
              (packet.artifactManifestEntries >= 5)
              (And
                (packet.profileEvidenceReady > 0)
                (packet.missingEvidence = 0)))))))

def contactParameterBenchValidationReadyNat
    (validation : ContactParameterBenchValidationNat) : Prop :=
  And
    (validation.packetReady > 0)
    (And
      (validation.fitRows > 0)
      (And
        (validation.holdoutRows > 0)
        (And
          (validation.completedMeasurements >= validation.fitRows + validation.holdoutRows)
          (And
            (validation.parameterResiduals > 0)
            (And
              (validation.fitPass > 0)
              (And
                (validation.holdoutPass > 0)
                (And
                  (validation.independentHoldoutPass > 0)
                  (validation.missingEvidence = 0))))))))

def contactParameterIntervalCalibrationReadyNat
    (calibration : ContactParameterIntervalCalibrationNat) : Prop :=
  And
    (calibration.benchValidationReady > 0)
    (And
      (calibration.parameterIntervals > 0)
      (And
        (calibration.acceptedParameterIntervals >= calibration.parameterIntervals)
        (And
          (calibration.uncertaintyRecords >= calibration.parameterIntervals)
          (And
            (calibration.holdoutAgreementRecords >= calibration.parameterIntervals)
            (And
              (calibration.simulatorParametersInsideBounds > 0)
              (calibration.missingEvidence = 0))))))

def externalPhysicsRunReadyNat
    (run : ExternalPhysicsRunNat) : Prop :=
  And
    (run.engineAvailable > 0)
    (And
      (run.modelExportReady > 0)
      (And
        (run.bodyRecords > 0)
        (And
          (run.resultBodyRecords >= run.bodyRecords)
          (And
            (run.stepsRequested > 0)
            (And
              (run.stepsCompleted >= run.stepsRequested)
              (run.missingEvidence = 0))))))

def externalPhysicsComparisonReadyNat
    (comparison : ExternalPhysicsComparisonNat) : Prop :=
  And
    (comparison.runReady > 0)
    (And
      (comparison.simulatorRecords > 0)
      (And
        (comparison.externalRecords >= comparison.simulatorRecords)
        (And
          (comparison.matchedRecords >= comparison.simulatorRecords)
          (And
            (comparison.toleranceRecords > 0)
            (And
              (comparison.withinTolerance > 0)
              (comparison.missingEvidence = 0))))))

def equilibriumRelationReadyNat
    (relation : EquilibriumRelationNat) : Prop :=
  And
    (relation.realizationReady > 0)
    (And
      (relation.solverSuccess > 0)
      (And
        (relation.contactStateReady > 0)
        (And
          (relation.energyTerms > 0)
          (And
            (relation.nonnegativeEnergyTerms >= relation.energyTerms)
            (And
              (relation.balanceResidual <= relation.residualTolerance)
              (relation.missingEvidence = 0))))))

def reachableEquilibriumControllabilityReadyNat
    (reachability : ReachableEquilibriumControllabilityNat) : Prop :=
  And
    (reachability.equilibriumReady > 0)
    (And
      (reachability.actuatorBasis > 0)
      (And
        (reachability.responseColumns > 0)
        (And
          (reachability.targetCells > 0)
          (And
            (reachability.reachableResponses > 0)
            (And
              (reachability.topologyEvidence > 0)
              (And
                (reachability.fullTargetRequired = 0 ∨
                  reachability.targetUnderactuated = 0)
                (reachability.missingEvidence = 0)))))))

def reachableEquilibriumBenchProtocolReadyNat
    (protocol : ReachableEquilibriumBenchProtocolNat) : Prop :=
  And
    (protocol.reachabilityReady > 0)
    (And
      (protocol.scenarioCount > 0)
      (And
        (protocol.actuatorColumnTrials > 0)
        (And
          (protocol.targetObservationCells > 0)
          (And
            (protocol.measurementColumns > 0)
            (And
              (protocol.topologyControlPolicy > 0)
              (And
                (protocol.passFailCriteria > 0)
                (protocol.missingEvidence = 0)))))))

def reachableEquilibriumBenchValidationReadyNat
    (validation : ReachableEquilibriumBenchValidationNat) : Prop :=
  And
    (validation.protocolReady > 0)
    (And
      (validation.resultRows > 0)
      (And
        (validation.targetMeasurements > 0)
        (And
          (validation.completedMeasurements >= validation.targetMeasurements)
          (And
            (validation.reachabilityChecks > 0)
            (And
              (validation.topologyLeakageChecks > 0)
              (And
                (validation.groupSequenceChecks > 0)
                (And
                  (validation.comparisonPass > 0)
                  (validation.missingEvidence = 0))))))))

def reachableEquilibriumAmplitudeCalibrationReadyNat
    (calibration : ReachableEquilibriumAmplitudeCalibrationNat) : Prop :=
  And
    (calibration.benchValidationReady > 0)
    (And
      (calibration.repeatedTrialGroups > 0)
      (And
        (calibration.amplitudeEstimates > 0)
        (And
          (calibration.residualFieldCells > 0)
          (And
            (calibration.uncertaintyBands > 0)
            (And
              (calibration.topologyLeakageBands > 0)
              (And
                (calibration.groupSequenceResiduals > 0)
                (calibration.missingEvidence = 0)))))))

def reachableEquilibriumEmpiricalProfileReadyNat
    (profile : ReachableEquilibriumEmpiricalProfileNat) : Prop :=
  And
    (profile.amplitudeCalibrationReady > 0)
    (And
      (profile.boundedProposals > 0)
      (And
        (profile.safeProposals > 0)
        (And
          (profile.toleranceRecords > 0)
          (And
            (profile.uncertaintyRecords > 0)
            (And
              (profile.holdoutHooks > 0)
              (profile.missingEvidence = 0))))))

def reachableEquilibriumProfileInverseReadyNat
    (inverse : ReachableEquilibriumProfileInverseNat) : Prop :=
  And
    (inverse.empiricalProfileReady > 0)
    (And
      (inverse.safeProfileProposals > 0)
      (And
        (inverse.targetCells > 0)
        (And
          (inverse.inverseSolveSuccess > 0)
          (And
            (inverse.residualRecords > 0)
            (And
              (inverse.scoreTerms > 0)
              (And
                (inverse.readOnlyProfileUse > 0)
                (inverse.missingEvidence = 0)))))))

def reachableEquilibriumProfileInverseAcceptanceReadyNat
    (acceptance : ReachableEquilibriumProfileInverseAcceptanceNat) : Prop :=
  And
    (acceptance.profileInverseReady > 0)
    (And
      (acceptance.scoreRecords > 0)
      (And
        (acceptance.residualScore <= acceptance.residualLimit)
        (And
          (acceptance.bandFailures <= acceptance.maxBandFailures)
          (And
            (acceptance.actuatorCount <= acceptance.actuatorLimit)
            (And
              (acceptance.acceptedForPreview > 0)
              (acceptance.missingEvidence = 0))))))

def reachableEquilibriumProfileInversePreviewPacketReadyNat
    (packet : ReachableEquilibriumProfileInversePreviewPacketNat) : Prop :=
  And
    (packet.acceptanceReady > 0)
    (And
      (packet.commandRecords > 0)
      (And
        (packet.eventRecords >= packet.commandRecords)
        (And
          (packet.targetRecords > 0)
          (And
            (packet.residualRecords > 0)
            (And
              (packet.readOnlyPacket > 0)
              (packet.missingEvidence = 0))))))

def reachableEquilibriumProfileInversePreviewReplayReadyNat
    (replay : ReachableEquilibriumProfileInversePreviewReplayNat) : Prop :=
  And
    (replay.packetReady > 0)
    (And
      (replay.commandRecords > 0)
      (And
        (replay.replayedCommands = replay.commandRecords)
        (And
          (replay.eventCoverage > 0)
          (And
            (replay.simulationRecords > 0)
            (And
              (replay.targetResidualRecords > 0)
              (And
                (replay.residualAgreement > 0)
                (And
                  (replay.readOnlyReplay > 0)
                  (replay.missingEvidence = 0))))))))

def reachableEquilibriumProfileInversePreviewPhysicalReadyNat
    (physical : ReachableEquilibriumProfileInversePreviewPhysicalNat) : Prop :=
  And
    (physical.replayReady > 0)
    (And
      (physical.physicalSolverSuccess > 0)
      (And
        (physical.commandRecords > 0)
        (And
          (physical.targetResidualRecords > 0)
          (And
            (physical.modelComparisonRecords > 0)
            (And
              (physical.energyRecords > 0)
              (And
                (physical.readOnlyPhysicalPreview > 0)
                (physical.missingEvidence = 0)))))))

def scaledSpringEnergyNat
    (stiffness displacement denominator : Nat) : ScaledNatQuantity :=
  { numerator := springEnergyNat stiffness displacement, denominator := denominator }

def scaledContactPenaltyNat
    (stiffness penetration denominator : Nat) : ScaledNatQuantity :=
  { numerator := contactPenaltyNat stiffness penetration, denominator := denominator }

def scaledSignedLoadWorkInt
    (load displacement : Int) (denominator : Nat) : ScaledIntQuantity :=
  { numerator := signedLoadWorkInt load displacement, denominator := denominator }

def scaledSignedWorkResidualInt
    (measured simulated : Int) (denominator : Nat) : ScaledNatQuantity :=
  { numerator := signedWorkResidualInt measured simulated, denominator := denominator }

def scaledSignedEnergyResidualTripleInt
    (measuredSigned simulatedSigned
      measuredMagnitude simulatedMagnitude
      measuredContact simulatedContact : Int)
    (denominator : Nat) : ScaledNatQuantity :=
  {
    numerator :=
      signedEnergyResidualTripleInt
        measuredSigned simulatedSigned
        measuredMagnitude simulatedMagnitude
        measuredContact simulatedContact,
    denominator := denominator,
  }

def measurementUnitScaleNat
    (scale : MeasurementUnitScale) (value : Nat) : ScaledNatQuantity :=
  { numerator := scale.numerator * value, denominator := scale.denominator }

def measurementUnitScaleInt
    (scale : MeasurementUnitScale) (value : Int) : ScaledIntQuantity :=
  { numerator := (scale.numerator : Int) * value, denominator := scale.denominator }

def measurementUnitScaleResidualNat
    (scale : MeasurementUnitScale) (measured simulated : Nat) :
    ScaledNatQuantity :=
  measurementUnitScaleNat scale (measuredWorkResidualNat measured simulated)

def measurementUnitScaleResidualInt
    (scale : MeasurementUnitScale) (measured simulated : Int) :
    ScaledNatQuantity :=
  measurementUnitScaleNat scale (signedWorkResidualInt measured simulated)

theorem springEnergyNat_nonnegative
    (stiffness displacement : Nat) :
    0 <= springEnergyNat stiffness displacement := by
  exact Nat.zero_le _

theorem hingeEnergyNat_nonnegative
    (stiffness angleError : Nat) :
    0 <= hingeEnergyNat stiffness angleError := by
  exact Nat.zero_le _

theorem contactPenaltyNat_nonnegative
    (stiffness penetration : Nat) :
    0 <= contactPenaltyNat stiffness penetration := by
  exact Nat.zero_le _

theorem scaledSpringEnergyNat_numerator_nonnegative
    (stiffness displacement denominator : Nat) :
    0 <= (scaledSpringEnergyNat stiffness displacement denominator).numerator := by
  exact Nat.zero_le _

theorem scaledSpringEnergyNat_zero_displacement
    (stiffness denominator : Nat) :
    (scaledSpringEnergyNat stiffness 0 denominator).numerator = 0 := by
  unfold scaledSpringEnergyNat springEnergyNat
  simp

theorem scaledSpringEnergyNat_preserves_denominator
    (stiffness displacement denominator : Nat) :
    (scaledSpringEnergyNat stiffness displacement denominator).denominator =
      denominator := by
  rfl

theorem scaledContactPenaltyNat_numerator_nonnegative
    (stiffness penetration denominator : Nat) :
    0 <= (scaledContactPenaltyNat stiffness penetration denominator).numerator := by
  exact Nat.zero_le _

theorem scaledContactPenaltyNat_zero_penetration
    (stiffness denominator : Nat) :
    (scaledContactPenaltyNat stiffness 0 denominator).numerator = 0 := by
  unfold scaledContactPenaltyNat contactPenaltyNat
  simp

theorem measurementUnitScaleNat_zero
    (scale : MeasurementUnitScale) :
    (measurementUnitScaleNat scale 0).numerator = 0 := by
  unfold measurementUnitScaleNat
  simp

theorem measurementUnitScaleNat_preserves_denominator
    (scale : MeasurementUnitScale) (value : Nat) :
    (measurementUnitScaleNat scale value).denominator = scale.denominator := by
  rfl

theorem clearanceExcessNat_zero_inside {clearance command : Nat}
    (hinside : command <= clearance) :
    clearanceExcessNat clearance command = 0 := by
  unfold clearanceExcessNat
  exact Nat.sub_eq_zero_of_le hinside

theorem verticalResidualStepNat_zero_inside {clearance gain command : Nat}
    (hinside : command <= clearance) :
    verticalResidualStepNat clearance gain command = 0 := by
  unfold verticalResidualStepNat
  rw [clearanceExcessNat_zero_inside hinside]
  simp

theorem contactPenaltyFromClearanceNat_nonnegative
    (stiffness clearance command : Nat) :
    0 <= contactPenaltyFromClearanceNat stiffness clearance command := by
  exact Nat.zero_le _

theorem loadWorkMagnitudeNat_nonnegative
    (load displacement : Nat) :
    0 <= loadWorkMagnitudeNat load displacement := by
  exact Nat.zero_le _

theorem loadWorkMagnitudeNat_zero_fixed
    (load : Nat) :
    loadWorkMagnitudeNat load 0 = 0 := by
  unfold loadWorkMagnitudeNat
  simp

theorem signedLoadWorkInt_zero_displacement
    (load : Int) :
    signedLoadWorkInt load 0 = 0 := by
  unfold signedLoadWorkInt
  simp

theorem signedLoadWorkInt_zero_load
    (displacement : Int) :
    signedLoadWorkInt 0 displacement = 0 := by
  unfold signedLoadWorkInt
  simp

theorem scaledSignedLoadWorkInt_zero_displacement
    (load : Int) (denominator : Nat) :
    (scaledSignedLoadWorkInt load 0 denominator).numerator = 0 := by
  unfold scaledSignedLoadWorkInt
  simp [signedLoadWorkInt_zero_displacement]

theorem scaledSignedLoadWorkInt_zero_load
    (displacement : Int) (denominator : Nat) :
    (scaledSignedLoadWorkInt 0 displacement denominator).numerator = 0 := by
  unfold scaledSignedLoadWorkInt
  simp [signedLoadWorkInt_zero_load]

theorem scaledSignedLoadWorkInt_preserves_denominator
    (load displacement : Int) (denominator : Nat) :
    (scaledSignedLoadWorkInt load displacement denominator).denominator =
      denominator := by
  rfl

theorem measurementUnitScaleInt_zero
    (scale : MeasurementUnitScale) :
    (measurementUnitScaleInt scale 0).numerator = 0 := by
  unfold measurementUnitScaleInt
  simp

theorem measurementUnitScaleInt_preserves_denominator
    (scale : MeasurementUnitScale) (value : Int) :
    (measurementUnitScaleInt scale value).denominator = scale.denominator := by
  rfl

theorem mechanicalStoredEnergyNat_nonnegative
    (spring hinge lock contact loadMagnitude : Nat) :
    0 <= mechanicalStoredEnergyNat spring hinge lock contact loadMagnitude := by
  exact Nat.zero_le _

theorem mechanicalStoredEnergyNat_zero_components :
    mechanicalStoredEnergyNat 0 0 0 0 0 = 0 := by
  unfold mechanicalStoredEnergyNat
  simp

theorem absoluteErrorNat_nonnegative
    (measured simulated : Nat) :
    0 <= absoluteErrorNat measured simulated := by
  exact Nat.zero_le _

theorem absoluteErrorNat_self
    (value : Nat) :
    absoluteErrorNat value value = 0 := by
  unfold absoluteErrorNat
  simp

theorem measuredWorkResidualNat_zero_when_equal
    (value : Nat) :
    measuredWorkResidualNat value value = 0 := by
  unfold measuredWorkResidualNat
  exact absoluteErrorNat_self value

theorem absoluteErrorInt_self
    (value : Int) :
    absoluteErrorInt value value = 0 := by
  unfold absoluteErrorInt
  simp

theorem signedWorkResidualInt_zero_when_equal
    (value : Int) :
    signedWorkResidualInt value value = 0 := by
  unfold signedWorkResidualInt
  exact absoluteErrorInt_self value

theorem scaledSignedWorkResidualInt_numerator_nonnegative
    (measured simulated : Int) (denominator : Nat) :
    0 <= (scaledSignedWorkResidualInt measured simulated denominator).numerator := by
  exact Nat.zero_le _

theorem scaledSignedWorkResidualInt_zero_when_equal
    (value : Int) (denominator : Nat) :
    (scaledSignedWorkResidualInt value value denominator).numerator = 0 := by
  unfold scaledSignedWorkResidualInt
  simp [signedWorkResidualInt_zero_when_equal]

theorem measurementUnitScaleResidualNat_zero_when_equal
    (scale : MeasurementUnitScale) (value : Nat) :
    (measurementUnitScaleResidualNat scale value value).numerator = 0 := by
  unfold measurementUnitScaleResidualNat
  simp [measuredWorkResidualNat_zero_when_equal, measurementUnitScaleNat_zero]

theorem measurementUnitScaleResidualNat_preserves_denominator
    (scale : MeasurementUnitScale) (measured simulated : Nat) :
    (measurementUnitScaleResidualNat scale measured simulated).denominator =
      scale.denominator := by
  unfold measurementUnitScaleResidualNat
  exact measurementUnitScaleNat_preserves_denominator scale
    (measuredWorkResidualNat measured simulated)

theorem measurementUnitScaleResidualInt_zero_when_equal
    (scale : MeasurementUnitScale) (value : Int) :
    (measurementUnitScaleResidualInt scale value value).numerator = 0 := by
  unfold measurementUnitScaleResidualInt
  simp [signedWorkResidualInt_zero_when_equal, measurementUnitScaleNat_zero]

theorem measurementUnitScaleResidualInt_preserves_denominator
    (scale : MeasurementUnitScale) (measured simulated : Int) :
    (measurementUnitScaleResidualInt scale measured simulated).denominator =
      scale.denominator := by
  unfold measurementUnitScaleResidualInt
  exact measurementUnitScaleNat_preserves_denominator scale
    (signedWorkResidualInt measured simulated)

theorem hardwareProfileCoverageComplete_true_when_equal
    (measured total : Nat) (heq : measured = total) :
    hardwareProfileCoverageComplete
      { measuredFields := measured, totalFields := total } = true := by
  unfold hardwareProfileCoverageComplete
  simp [heq]

theorem hardwareProfileCoverageMissing_zero_when_complete
    (measured total : Nat) (heq : measured = total) :
    hardwareProfileCoverageMissing
      { measuredFields := measured, totalFields := total } = 0 := by
  unfold hardwareProfileCoverageMissing
  simp [heq]

theorem hardwareProfileCoverageMissing_preserves_total
    (measured total : Nat) :
    (hardwareProfileCoverageMissing
      { measuredFields := measured, totalFields := total }) =
      total - measured := by
  rfl

theorem calibrationFitResidualPassNat_zero
    {sampleCount tolerance : Nat} (hcount : sampleCount > 0) :
    calibrationFitResidualPassNat
      {
        sampleCount := sampleCount,
        rawResidual := 0,
        fittedResidual := 0,
        tolerance := tolerance,
      } := by
  unfold calibrationFitResidualPassNat
  constructor
  · exact hcount
  · exact Nat.zero_le tolerance

theorem calibrationFitResidualNat_preserves_sample_count
    (fit : CalibrationFitResidualNat) :
    fit.sampleCount = fit.sampleCount := by
  rfl

theorem calibrationModelProfileUpdateSafeNat_intro
    (update : CalibrationModelProfileUpdateNat)
    (hcount : update.sampleCount > 0)
    (hlower : update.lowerBound <= update.proposed)
    (hupper : update.proposed <= update.upperBound) :
    calibrationModelProfileUpdateSafeNat update := by
  unfold calibrationModelProfileUpdateSafeNat
  exact And.intro hcount (And.intro hlower hupper)

theorem calibrationModelProfileUpdateSafeNat_bounds
    {update : CalibrationModelProfileUpdateNat}
    (hsafe : calibrationModelProfileUpdateSafeNat update) :
    update.lowerBound <= update.proposed ∧
      update.proposed <= update.upperBound := by
  exact hsafe.2

theorem calibrationModelProfileCandidateImprovesNat_intro
    (candidate : CalibrationModelProfileCandidateNat)
    (happlied : candidate.appliedUpdates > 0)
    (hresidual : candidate.afterResidual <= candidate.beforeResidual)
    (hmissing : candidate.missingAfter <= candidate.missingBefore) :
    calibrationModelProfileCandidateImprovesNat candidate := by
  unfold calibrationModelProfileCandidateImprovesNat
  exact And.intro happlied (And.intro hresidual hmissing)

theorem calibrationModelProfileCandidateImprovesNat_applied_updates
    {candidate : CalibrationModelProfileCandidateNat}
    (himproves : calibrationModelProfileCandidateImprovesNat candidate) :
    candidate.appliedUpdates > 0 := by
  exact himproves.1

theorem calibrationModelProfileCandidateImprovesNat_residual_nonincrease
    {candidate : CalibrationModelProfileCandidateNat}
    (himproves : calibrationModelProfileCandidateImprovesNat candidate) :
    candidate.afterResidual <= candidate.beforeResidual := by
  exact himproves.2.1

theorem calibrationModelProfileCandidateImprovesNat_missing_nonincrease
    {candidate : CalibrationModelProfileCandidateNat}
    (himproves : calibrationModelProfileCandidateImprovesNat candidate) :
    candidate.missingAfter <= candidate.missingBefore := by
  exact himproves.2.2

theorem calibrationModelProfileCandidateScoreNat_le_before
    {candidate : CalibrationModelProfileCandidateNat}
    (himproves : calibrationModelProfileCandidateImprovesNat candidate) :
    calibrationModelProfileCandidateScoreNat candidate <=
      candidate.beforeResidual := by
  unfold calibrationModelProfileCandidateScoreNat
  exact himproves.2.1

theorem calibrationModelProfileHoldoutPassNat_intro
    (validation : CalibrationModelProfileHoldoutNat)
    (happlied : validation.appliedUpdates > 0)
    (hfit : validation.fitAfterResidual <= validation.fitBeforeResidual)
    (hholdout : validation.holdoutAfterResidual <= validation.holdoutBeforeResidual)
    (hfitMissing : validation.fitMissingAfter <= validation.fitMissingBefore)
    (hholdoutMissing : validation.holdoutMissingAfter <= validation.holdoutMissingBefore) :
    calibrationModelProfileHoldoutPassNat validation := by
  unfold calibrationModelProfileHoldoutPassNat
  exact And.intro happlied
    (And.intro hfit
      (And.intro hholdout
        (And.intro hfitMissing hholdoutMissing)))

theorem calibrationModelProfileHoldoutPassNat_applied_updates
    {validation : CalibrationModelProfileHoldoutNat}
    (hpass : calibrationModelProfileHoldoutPassNat validation) :
    validation.appliedUpdates > 0 := by
  exact hpass.1

theorem calibrationModelProfileHoldoutPassNat_fit_residual_nonincrease
    {validation : CalibrationModelProfileHoldoutNat}
    (hpass : calibrationModelProfileHoldoutPassNat validation) :
    validation.fitAfterResidual <= validation.fitBeforeResidual := by
  exact hpass.2.1

theorem calibrationModelProfileHoldoutPassNat_holdout_residual_nonincrease
    {validation : CalibrationModelProfileHoldoutNat}
    (hpass : calibrationModelProfileHoldoutPassNat validation) :
    validation.holdoutAfterResidual <= validation.holdoutBeforeResidual := by
  exact hpass.2.2.1

theorem calibrationModelProfileHoldoutPassNat_missing_nonincrease
    {validation : CalibrationModelProfileHoldoutNat}
    (hpass : calibrationModelProfileHoldoutPassNat validation) :
    validation.fitMissingAfter <= validation.fitMissingBefore ∧
      validation.holdoutMissingAfter <= validation.holdoutMissingBefore := by
  exact hpass.2.2.2

theorem calibrationModelProfileHoldoutScoreNat_le_before
    {validation : CalibrationModelProfileHoldoutNat}
    (hpass : calibrationModelProfileHoldoutPassNat validation) :
    calibrationModelProfileHoldoutScoreNat validation <=
      validation.holdoutBeforeResidual := by
  unfold calibrationModelProfileHoldoutScoreNat
  exact hpass.2.2.1

theorem calibrationTrainHoldoutSplitReadyNat_intro
    (split : CalibrationTrainHoldoutSplitNat)
    (hfit : split.fitSamples > 0)
    (hholdout : split.holdoutSamples > 0)
    (hoverlap : split.overlapSamples = 0)
    (hfrozen : split.frozenProfile > 0) :
    calibrationTrainHoldoutSplitReadyNat split := by
  unfold calibrationTrainHoldoutSplitReadyNat
  exact And.intro hfit (And.intro hholdout (And.intro hoverlap hfrozen))

theorem calibrationTrainHoldoutSplitReadyNat_has_fit_samples
    {split : CalibrationTrainHoldoutSplitNat}
    (hready : calibrationTrainHoldoutSplitReadyNat split) :
    split.fitSamples > 0 := by
  exact hready.1

theorem calibrationTrainHoldoutSplitReadyNat_has_holdout_samples
    {split : CalibrationTrainHoldoutSplitNat}
    (hready : calibrationTrainHoldoutSplitReadyNat split) :
    split.holdoutSamples > 0 := by
  exact hready.2.1

theorem calibrationTrainHoldoutSplitReadyNat_zero_overlap
    {split : CalibrationTrainHoldoutSplitNat}
    (hready : calibrationTrainHoldoutSplitReadyNat split) :
    split.overlapSamples = 0 := by
  exact hready.2.2.1

theorem calibrationTrainHoldoutSplitReadyNat_profile_frozen
    {split : CalibrationTrainHoldoutSplitNat}
    (hready : calibrationTrainHoldoutSplitReadyNat split) :
    split.frozenProfile > 0 := by
  exact hready.2.2.2

theorem calibrationTrainHoldoutProvenanceReadyNat_intro
    (provenance : CalibrationTrainHoldoutProvenanceNat)
    (hfitId : provenance.fitDatasetIdKnown > 0)
    (hholdoutId : provenance.holdoutDatasetIdKnown > 0)
    (hdatasets : provenance.distinctDatasetIds > 0)
    (hfiles : provenance.distinctSourceFileIds > 0)
    (hfitRole : provenance.fitRoleMarked > 0)
    (hholdoutRole : provenance.holdoutRoleMarked > 0)
    (hfrozen : provenance.profileFrozen > 0)
    (hprofile : provenance.holdoutProfileMatchesFrozenProfile > 0) :
    calibrationTrainHoldoutProvenanceReadyNat provenance := by
  unfold calibrationTrainHoldoutProvenanceReadyNat
  exact And.intro hfitId
    (And.intro hholdoutId
      (And.intro hdatasets
        (And.intro hfiles
          (And.intro hfitRole
            (And.intro hholdoutRole
              (And.intro hfrozen hprofile))))))

theorem calibrationTrainHoldoutProvenanceReadyNat_distinct_files
    {provenance : CalibrationTrainHoldoutProvenanceNat}
    (hready : calibrationTrainHoldoutProvenanceReadyNat provenance) :
    provenance.distinctSourceFileIds > 0 := by
  exact hready.2.2.2.1

theorem calibrationTrainHoldoutProvenanceReadyNat_profile_frozen
    {provenance : CalibrationTrainHoldoutProvenanceNat}
    (hready : calibrationTrainHoldoutProvenanceReadyNat provenance) :
    provenance.profileFrozen > 0 := by
  exact hready.2.2.2.2.2.2.1

theorem calibrationTrainHoldoutProvenanceReadyNat_profile_matches
    {provenance : CalibrationTrainHoldoutProvenanceNat}
    (hready : calibrationTrainHoldoutProvenanceReadyNat provenance) :
    provenance.holdoutProfileMatchesFrozenProfile > 0 := by
  exact hready.2.2.2.2.2.2.2

theorem calibrationBenchProtocolCoverageReadyNat_intro
    (coverage : CalibrationBenchProtocolCoverageNat)
    (hscenario : coverage.scenarioCount > 0)
    (hinstrument : coverage.instrumentCount > 0)
    (houtput : coverage.outputArtifactCount > 0)
    (hcolumns : coverage.measurementColumnCount > 0)
    (hroles : coverage.datasetRoleCount >= 2) :
    calibrationBenchProtocolCoverageReadyNat coverage := by
  unfold calibrationBenchProtocolCoverageReadyNat
  exact And.intro hscenario
    (And.intro hinstrument
      (And.intro houtput
        (And.intro hcolumns hroles)))

theorem calibrationBenchProtocolCoverageReadyNat_has_scenarios
    {coverage : CalibrationBenchProtocolCoverageNat}
    (hready : calibrationBenchProtocolCoverageReadyNat coverage) :
    coverage.scenarioCount > 0 := by
  exact hready.1

theorem calibrationBenchProtocolCoverageReadyNat_has_outputs
    {coverage : CalibrationBenchProtocolCoverageNat}
    (hready : calibrationBenchProtocolCoverageReadyNat coverage) :
    coverage.outputArtifactCount > 0 := by
  exact hready.2.2.1

theorem calibrationBenchProtocolCoverageReadyNat_has_two_dataset_roles
    {coverage : CalibrationBenchProtocolCoverageNat}
    (hready : calibrationBenchProtocolCoverageReadyNat coverage) :
    coverage.datasetRoleCount >= 2 := by
  exact hready.2.2.2.2

theorem calibrationBenchPacketCompleteNat_intro
    (packet : CalibrationBenchPacketCompletenessNat)
    (hnotebook : packet.notebookArtifacts > 0)
    (hprotocol : packet.protocolArtifacts > 0)
    (hfit : packet.fitTemplateArtifacts > 0)
    (hholdout : packet.holdoutTemplateArtifacts > 0)
    (hvalidation : packet.validationInstructionArtifacts > 0)
    (hmanifest : packet.manifestEntries >= 5) :
    calibrationBenchPacketCompleteNat packet := by
  unfold calibrationBenchPacketCompleteNat
  exact And.intro hnotebook
    (And.intro hprotocol
      (And.intro hfit
        (And.intro hholdout
          (And.intro hvalidation hmanifest))))

theorem calibrationBenchPacketCompleteNat_has_notebook
    {packet : CalibrationBenchPacketCompletenessNat}
    (hcomplete : calibrationBenchPacketCompleteNat packet) :
    packet.notebookArtifacts > 0 := by
  exact hcomplete.1

theorem calibrationBenchPacketCompleteNat_has_fit_template
    {packet : CalibrationBenchPacketCompletenessNat}
    (hcomplete : calibrationBenchPacketCompleteNat packet) :
    packet.fitTemplateArtifacts > 0 := by
  exact hcomplete.2.2.1

theorem calibrationBenchPacketCompleteNat_has_holdout_template
    {packet : CalibrationBenchPacketCompletenessNat}
    (hcomplete : calibrationBenchPacketCompleteNat packet) :
    packet.holdoutTemplateArtifacts > 0 := by
  exact hcomplete.2.2.2.1

theorem calibrationBenchPacketCompleteNat_has_manifest
    {packet : CalibrationBenchPacketCompletenessNat}
    (hcomplete : calibrationBenchPacketCompleteNat packet) :
    packet.manifestEntries >= 5 := by
  exact hcomplete.2.2.2.2.2

theorem calibrationBenchExecutedValidationReadyNat_intro
    (validation : CalibrationBenchExecutedValidationNat)
    (hfit : validation.fitSamples > 0)
    (hholdout : validation.holdoutSamples > 0)
    (hupdates : validation.appliedUpdates > 0)
    (hresidual : validation.residualValidationPass > 0)
    (hindependent : validation.independentValidationPass > 0)
    (hmissing : validation.missingEvidence = 0) :
    calibrationBenchExecutedValidationReadyNat validation := by
  unfold calibrationBenchExecutedValidationReadyNat
  exact And.intro hfit
    (And.intro hholdout
      (And.intro hupdates
        (And.intro hresidual
          (And.intro hindependent hmissing))))

theorem calibrationBenchExecutedValidationReadyNat_has_fit_samples
    {validation : CalibrationBenchExecutedValidationNat}
    (hready : calibrationBenchExecutedValidationReadyNat validation) :
    validation.fitSamples > 0 := by
  exact hready.1

theorem calibrationBenchExecutedValidationReadyNat_has_holdout_samples
    {validation : CalibrationBenchExecutedValidationNat}
    (hready : calibrationBenchExecutedValidationReadyNat validation) :
    validation.holdoutSamples > 0 := by
  exact hready.2.1

theorem calibrationBenchExecutedValidationReadyNat_has_applied_updates
    {validation : CalibrationBenchExecutedValidationNat}
    (hready : calibrationBenchExecutedValidationReadyNat validation) :
    validation.appliedUpdates > 0 := by
  exact hready.2.2.1

theorem calibrationBenchExecutedValidationReadyNat_has_independent_validation
    {validation : CalibrationBenchExecutedValidationNat}
    (hready : calibrationBenchExecutedValidationReadyNat validation) :
    validation.independentValidationPass > 0 := by
  exact hready.2.2.2.2.1

theorem calibrationBenchExecutedValidationReadyNat_zero_missing_evidence
    {validation : CalibrationBenchExecutedValidationNat}
    (hready : calibrationBenchExecutedValidationReadyNat validation) :
    validation.missingEvidence = 0 := by
  exact hready.2.2.2.2.2

theorem physicalValidationReadyNat_intro
    (readiness : PhysicalValidationReadinessNat)
    (hcalibration : readiness.calibrationExecutionReady > 0)
    (hscenarios : readiness.verticalLoadScenarios > 0)
    (hvertical : readiness.verticalLoadComparisonPass > 0)
    (hload : readiness.loadProxyTerms > 0)
    (hcontact : readiness.contactProxyTerms > 0)
    (hclearance : readiness.clearanceConfigured > 0)
    (hmissing : readiness.missingEvidence = 0) :
    physicalValidationReadyNat readiness := by
  unfold physicalValidationReadyNat
  exact And.intro hcalibration
    (And.intro hscenarios
      (And.intro hvertical
        (And.intro hload
          (And.intro hcontact
            (And.intro hclearance hmissing)))))

theorem physicalValidationReadyNat_has_calibration
    {readiness : PhysicalValidationReadinessNat}
    (hready : physicalValidationReadyNat readiness) :
    readiness.calibrationExecutionReady > 0 := by
  exact hready.1

theorem physicalValidationReadyNat_has_vertical_load
    {readiness : PhysicalValidationReadinessNat}
    (hready : physicalValidationReadyNat readiness) :
    readiness.verticalLoadComparisonPass > 0 := by
  exact hready.2.2.1

theorem physicalValidationReadyNat_has_load_proxy
    {readiness : PhysicalValidationReadinessNat}
    (hready : physicalValidationReadyNat readiness) :
    readiness.loadProxyTerms > 0 := by
  exact hready.2.2.2.1

theorem physicalValidationReadyNat_has_contact_proxy
    {readiness : PhysicalValidationReadinessNat}
    (hready : physicalValidationReadyNat readiness) :
    readiness.contactProxyTerms > 0 := by
  exact hready.2.2.2.2.1

theorem physicalValidationReadyNat_zero_missing_evidence
    {readiness : PhysicalValidationReadinessNat}
    (hready : physicalValidationReadyNat readiness) :
    readiness.missingEvidence = 0 := by
  exact hready.2.2.2.2.2.2

theorem contactStateAbstractionReadyNat_intro
    (contact : ContactStateAbstractionNat)
    (hbodies : contact.bodyCount > 0)
    (hpins : contact.pinCount = contact.bodyCount)
    (hholes : contact.holeCount = contact.bodyCount)
    (hclearancePairs : contact.clearancePairCount = contact.bodyCount)
    (hrecords : contact.contactStateRecords >= contact.bodyCount)
    (hpenalty : contact.penaltyTerms >= contact.bodyCount)
    (hclearance : contact.clearanceConfigured > 0)
    (hmissing : contact.missingEvidence = 0) :
    contactStateAbstractionReadyNat contact := by
  unfold contactStateAbstractionReadyNat
  exact And.intro hbodies
    (And.intro hpins
      (And.intro hholes
        (And.intro hclearancePairs
          (And.intro hrecords
            (And.intro hpenalty
              (And.intro hclearance hmissing))))))

theorem contactStateAbstractionReadyNat_has_bodies
    {contact : ContactStateAbstractionNat}
    (hready : contactStateAbstractionReadyNat contact) :
    contact.bodyCount > 0 := by
  exact hready.1

theorem contactStateAbstractionReadyNat_pins_match_bodies
    {contact : ContactStateAbstractionNat}
    (hready : contactStateAbstractionReadyNat contact) :
    contact.pinCount = contact.bodyCount := by
  exact hready.2.1

theorem contactStateAbstractionReadyNat_holes_match_bodies
    {contact : ContactStateAbstractionNat}
    (hready : contactStateAbstractionReadyNat contact) :
    contact.holeCount = contact.bodyCount := by
  exact hready.2.2.1

theorem contactStateAbstractionReadyNat_has_contact_records
    {contact : ContactStateAbstractionNat}
    (hready : contactStateAbstractionReadyNat contact) :
    contact.contactStateRecords >= contact.bodyCount := by
  exact hready.2.2.2.2.1

theorem contactStateAbstractionReadyNat_has_penalty_terms
    {contact : ContactStateAbstractionNat}
    (hready : contactStateAbstractionReadyNat contact) :
    contact.penaltyTerms >= contact.bodyCount := by
  exact hready.2.2.2.2.2.1

theorem contactStateAbstractionReadyNat_zero_missing_evidence
    {contact : ContactStateAbstractionNat}
    (hready : contactStateAbstractionReadyNat contact) :
    contact.missingEvidence = 0 := by
  exact hready.2.2.2.2.2.2.2

theorem contactGraphConsistentNat_intro
    (graph : ContactGraphConsistencyNat)
    (hbodies : graph.activeBodies > 0)
    (hrecords : graph.contactRecords >= graph.activeBodies)
    (hedges : graph.removedIncidentActiveEdges = 0)
    (hcontacts : graph.removedActiveContacts = 0)
    (hsupport : graph.groupSupportRecords >= graph.groupSupportCells)
    (hmissing : graph.missingEvidence = 0) :
    contactGraphConsistentNat graph := by
  unfold contactGraphConsistentNat
  exact And.intro hbodies
    (And.intro hrecords
      (And.intro hedges
        (And.intro hcontacts
          (And.intro hsupport hmissing))))

theorem contactGraphConsistentNat_has_active_bodies
    {graph : ContactGraphConsistencyNat}
    (hready : contactGraphConsistentNat graph) :
    graph.activeBodies > 0 := by
  exact hready.1

theorem contactGraphConsistentNat_contact_records_cover_bodies
    {graph : ContactGraphConsistencyNat}
    (hready : contactGraphConsistentNat graph) :
    graph.contactRecords >= graph.activeBodies := by
  exact hready.2.1

theorem contactGraphConsistentNat_removed_edges_clear
    {graph : ContactGraphConsistencyNat}
    (hready : contactGraphConsistentNat graph) :
    graph.removedIncidentActiveEdges = 0 := by
  exact hready.2.2.1

theorem contactGraphConsistentNat_removed_contacts_clear
    {graph : ContactGraphConsistencyNat}
    (hready : contactGraphConsistentNat graph) :
    graph.removedActiveContacts = 0 := by
  exact hready.2.2.2.1

theorem contactGraphConsistentNat_support_records_cover_support
    {graph : ContactGraphConsistencyNat}
    (hready : contactGraphConsistentNat graph) :
    graph.groupSupportRecords >= graph.groupSupportCells := by
  exact hready.2.2.2.2.1

theorem contactGraphConsistentNat_zero_missing_evidence
    {graph : ContactGraphConsistencyNat}
    (hready : contactGraphConsistentNat graph) :
    graph.missingEvidence = 0 := by
  exact hready.2.2.2.2.2

theorem physicalRealizationMapReadyNat_intro
    (realization : PhysicalRealizationMapNat)
    (hoperators : realization.abstractOperators > 0)
    (hrealized : realization.realizedOperators >= realization.abstractOperators)
    (hsupport : realization.supportRecords >= realization.abstractOperators)
    (heffects : realization.stateEffectRecords >= realization.abstractOperators)
    (hlabels : realization.claimLabels >= realization.abstractOperators)
    (hcontact : realization.contactGraphReady > 0)
    (hmissing : realization.missingEvidence = 0) :
    physicalRealizationMapReadyNat realization := by
  unfold physicalRealizationMapReadyNat
  exact And.intro hoperators
    (And.intro hrealized
      (And.intro hsupport
        (And.intro heffects
          (And.intro hlabels
            (And.intro hcontact hmissing)))))

theorem physicalRealizationMapReadyNat_has_operators
    {realization : PhysicalRealizationMapNat}
    (hready : physicalRealizationMapReadyNat realization) :
    realization.abstractOperators > 0 := by
  exact hready.1

theorem physicalRealizationMapReadyNat_realizes_all
    {realization : PhysicalRealizationMapNat}
    (hready : physicalRealizationMapReadyNat realization) :
    realization.realizedOperators >= realization.abstractOperators := by
  exact hready.2.1

theorem physicalRealizationMapReadyNat_has_support_records
    {realization : PhysicalRealizationMapNat}
    (hready : physicalRealizationMapReadyNat realization) :
    realization.supportRecords >= realization.abstractOperators := by
  exact hready.2.2.1

theorem physicalRealizationMapReadyNat_has_state_effects
    {realization : PhysicalRealizationMapNat}
    (hready : physicalRealizationMapReadyNat realization) :
    realization.stateEffectRecords >= realization.abstractOperators := by
  exact hready.2.2.2.1

theorem physicalRealizationMapReadyNat_has_claim_labels
    {realization : PhysicalRealizationMapNat}
    (hready : physicalRealizationMapReadyNat realization) :
    realization.claimLabels >= realization.abstractOperators := by
  exact hready.2.2.2.2.1

theorem physicalRealizationMapReadyNat_has_contact_graph
    {realization : PhysicalRealizationMapNat}
    (hready : physicalRealizationMapReadyNat realization) :
    realization.contactGraphReady > 0 := by
  exact hready.2.2.2.2.2.1

theorem physicalRealizationMapReadyNat_zero_missing_evidence
    {realization : PhysicalRealizationMapNat}
    (hready : physicalRealizationMapReadyNat realization) :
    realization.missingEvidence = 0 := by
  exact hready.2.2.2.2.2.2

theorem externalPhysicsEngineAuditReadyNat_intro
    (audit : ExternalPhysicsEngineAuditNat)
    (hengines : audit.engineCandidates > 0)
    (havailable : audit.availableCandidates > 0)
    (hfeatures : audit.requiredFeatureRecords > 0)
    (hscenarios : audit.scenarioRecords > 0)
    (hcontacts : audit.contactModelRecords > 0)
    (hindependent : audit.independentToolRecords > 0)
    (hmissing : audit.missingEvidence = 0) :
    externalPhysicsEngineAuditReadyNat audit := by
  unfold externalPhysicsEngineAuditReadyNat
  exact And.intro hengines
    (And.intro havailable
      (And.intro hfeatures
        (And.intro hscenarios
          (And.intro hcontacts
            (And.intro hindependent hmissing)))))

theorem externalPhysicsEngineAuditReadyNat_has_engines
    {audit : ExternalPhysicsEngineAuditNat}
    (hready : externalPhysicsEngineAuditReadyNat audit) :
    audit.engineCandidates > 0 := by
  exact hready.1

theorem externalPhysicsEngineAuditReadyNat_has_available
    {audit : ExternalPhysicsEngineAuditNat}
    (hready : externalPhysicsEngineAuditReadyNat audit) :
    audit.availableCandidates > 0 := by
  exact hready.2.1

theorem externalPhysicsEngineAuditReadyNat_has_features
    {audit : ExternalPhysicsEngineAuditNat}
    (hready : externalPhysicsEngineAuditReadyNat audit) :
    audit.requiredFeatureRecords > 0 := by
  exact hready.2.2.1

theorem externalPhysicsEngineAuditReadyNat_has_scenarios
    {audit : ExternalPhysicsEngineAuditNat}
    (hready : externalPhysicsEngineAuditReadyNat audit) :
    audit.scenarioRecords > 0 := by
  exact hready.2.2.2.1

theorem externalPhysicsEngineAuditReadyNat_has_contact_model
    {audit : ExternalPhysicsEngineAuditNat}
    (hready : externalPhysicsEngineAuditReadyNat audit) :
    audit.contactModelRecords > 0 := by
  exact hready.2.2.2.2.1

theorem externalPhysicsEngineAuditReadyNat_has_independent_tool
    {audit : ExternalPhysicsEngineAuditNat}
    (hready : externalPhysicsEngineAuditReadyNat audit) :
    audit.independentToolRecords > 0 := by
  exact hready.2.2.2.2.2.1

theorem externalPhysicsEngineAuditReadyNat_zero_missing_evidence
    {audit : ExternalPhysicsEngineAuditNat}
    (hready : externalPhysicsEngineAuditReadyNat audit) :
    audit.missingEvidence = 0 := by
  exact hready.2.2.2.2.2.2

theorem externalPhysicsModelExportReadyNat_intro
    (modelExport : ExternalPhysicsModelExportNat)
    (hmodel : modelExport.modelRecords > 0)
    (hbodies : modelExport.bodyRecords > 0)
    (hfixed : modelExport.fixedBodyRecords > 0)
    (hgravity : modelExport.gravityRecords > 0)
    (hxml : modelExport.xmlBytes > 0)
    (hmissing : modelExport.missingEvidence = 0) :
    externalPhysicsModelExportReadyNat modelExport := by
  unfold externalPhysicsModelExportReadyNat
  exact And.intro hmodel
    (And.intro hbodies
      (And.intro hfixed
        (And.intro hgravity
          (And.intro hxml hmissing))))

theorem externalPhysicsModelExportReadyNat_has_model
    {modelExport : ExternalPhysicsModelExportNat}
    (hready : externalPhysicsModelExportReadyNat modelExport) :
    modelExport.modelRecords > 0 := by
  exact hready.1

theorem externalPhysicsModelExportReadyNat_has_bodies
    {modelExport : ExternalPhysicsModelExportNat}
    (hready : externalPhysicsModelExportReadyNat modelExport) :
    modelExport.bodyRecords > 0 := by
  exact hready.2.1

theorem externalPhysicsModelExportReadyNat_has_fixed_bodies
    {modelExport : ExternalPhysicsModelExportNat}
    (hready : externalPhysicsModelExportReadyNat modelExport) :
    modelExport.fixedBodyRecords > 0 := by
  exact hready.2.2.1

theorem externalPhysicsModelExportReadyNat_has_gravity
    {modelExport : ExternalPhysicsModelExportNat}
    (hready : externalPhysicsModelExportReadyNat modelExport) :
    modelExport.gravityRecords > 0 := by
  exact hready.2.2.2.1

theorem externalPhysicsModelExportReadyNat_has_xml
    {modelExport : ExternalPhysicsModelExportNat}
    (hready : externalPhysicsModelExportReadyNat modelExport) :
    modelExport.xmlBytes > 0 := by
  exact hready.2.2.2.2.1

theorem externalPhysicsModelExportReadyNat_zero_missing_evidence
    {modelExport : ExternalPhysicsModelExportNat}
    (hready : externalPhysicsModelExportReadyNat modelExport) :
    modelExport.missingEvidence = 0 := by
  exact hready.2.2.2.2.2

theorem externalContactGeometryReadyNat_intro
    (geometry : ExternalContactGeometryNat)
    (hpins : geometry.pinRecords > 0)
    (hholes : geometry.holeRecords >= geometry.pinRecords)
    (hclearance : geometry.clearanceRecords >= geometry.pinRecords)
    (hpairs : geometry.contactPairRecords >= geometry.pinRecords)
    (hactive : geometry.activeContactPairs >= geometry.pinRecords)
    (hxml : geometry.xmlBytes > 0)
    (hmissing : geometry.missingEvidence = 0) :
    externalContactGeometryReadyNat geometry := by
  unfold externalContactGeometryReadyNat
  exact And.intro hpins
    (And.intro hholes
      (And.intro hclearance
        (And.intro hpairs
          (And.intro hactive
            (And.intro hxml hmissing)))))

theorem externalContactGeometryReadyNat_has_pins
    {geometry : ExternalContactGeometryNat}
    (hready : externalContactGeometryReadyNat geometry) :
    geometry.pinRecords > 0 := by
  exact hready.1

theorem externalContactGeometryReadyNat_holes_cover_pins
    {geometry : ExternalContactGeometryNat}
    (hready : externalContactGeometryReadyNat geometry) :
    geometry.holeRecords >= geometry.pinRecords := by
  exact hready.2.1

theorem externalContactGeometryReadyNat_clearance_covers_pins
    {geometry : ExternalContactGeometryNat}
    (hready : externalContactGeometryReadyNat geometry) :
    geometry.clearanceRecords >= geometry.pinRecords := by
  exact hready.2.2.1

theorem externalContactGeometryReadyNat_pairs_cover_pins
    {geometry : ExternalContactGeometryNat}
    (hready : externalContactGeometryReadyNat geometry) :
    geometry.contactPairRecords >= geometry.pinRecords := by
  exact hready.2.2.2.1

theorem externalContactGeometryReadyNat_active_pairs_cover_pins
    {geometry : ExternalContactGeometryNat}
    (hready : externalContactGeometryReadyNat geometry) :
    geometry.activeContactPairs >= geometry.pinRecords := by
  exact hready.2.2.2.2.1

theorem externalContactGeometryReadyNat_has_xml
    {geometry : ExternalContactGeometryNat}
    (hready : externalContactGeometryReadyNat geometry) :
    geometry.xmlBytes > 0 := by
  exact hready.2.2.2.2.2.1

theorem externalContactGeometryReadyNat_zero_missing_evidence
    {geometry : ExternalContactGeometryNat}
    (hready : externalContactGeometryReadyNat geometry) :
    geometry.missingEvidence = 0 := by
  exact hready.2.2.2.2.2.2

theorem externalContactParameterProfileReadyNat_intro
    (profile : ExternalContactParameterProfileNat)
    (hpairs : profile.contactPairRecords > 0)
    (hparameters : profile.parameterRecords >= profile.contactPairRecords)
    (hfriction : profile.frictionRecords >= profile.contactPairRecords)
    (hsolver : profile.solverParameterRecords >= profile.contactPairRecords)
    (hstiffness : profile.stiffnessRecords >= profile.contactPairRecords)
    (hdamping : profile.dampingRecords >= profile.contactPairRecords)
    (hxml : profile.xmlAttributeBytes > 0)
    (hmissing : profile.missingEvidence = 0) :
    externalContactParameterProfileReadyNat profile := by
  unfold externalContactParameterProfileReadyNat
  exact And.intro hpairs
    (And.intro hparameters
      (And.intro hfriction
        (And.intro hsolver
          (And.intro hstiffness
            (And.intro hdamping
              (And.intro hxml hmissing))))))

theorem externalContactParameterProfileReadyNat_has_pairs
    {profile : ExternalContactParameterProfileNat}
    (hready : externalContactParameterProfileReadyNat profile) :
    profile.contactPairRecords > 0 := by
  exact hready.1

theorem externalContactParameterProfileReadyNat_parameters_cover_pairs
    {profile : ExternalContactParameterProfileNat}
    (hready : externalContactParameterProfileReadyNat profile) :
    profile.parameterRecords >= profile.contactPairRecords := by
  exact hready.2.1

theorem externalContactParameterProfileReadyNat_friction_covers_pairs
    {profile : ExternalContactParameterProfileNat}
    (hready : externalContactParameterProfileReadyNat profile) :
    profile.frictionRecords >= profile.contactPairRecords := by
  exact hready.2.2.1

theorem externalContactParameterProfileReadyNat_solver_covers_pairs
    {profile : ExternalContactParameterProfileNat}
    (hready : externalContactParameterProfileReadyNat profile) :
    profile.solverParameterRecords >= profile.contactPairRecords := by
  exact hready.2.2.2.1

theorem externalContactParameterProfileReadyNat_stiffness_covers_pairs
    {profile : ExternalContactParameterProfileNat}
    (hready : externalContactParameterProfileReadyNat profile) :
    profile.stiffnessRecords >= profile.contactPairRecords := by
  exact hready.2.2.2.2.1

theorem externalContactParameterProfileReadyNat_damping_covers_pairs
    {profile : ExternalContactParameterProfileNat}
    (hready : externalContactParameterProfileReadyNat profile) :
    profile.dampingRecords >= profile.contactPairRecords := by
  exact hready.2.2.2.2.2.1

theorem externalContactParameterProfileReadyNat_has_xml_attributes
    {profile : ExternalContactParameterProfileNat}
    (hready : externalContactParameterProfileReadyNat profile) :
    profile.xmlAttributeBytes > 0 := by
  exact hready.2.2.2.2.2.2.1

theorem externalContactParameterProfileReadyNat_zero_missing_evidence
    {profile : ExternalContactParameterProfileNat}
    (hready : externalContactParameterProfileReadyNat profile) :
    profile.missingEvidence = 0 := by
  exact hready.2.2.2.2.2.2.2

theorem contactParameterCalibrationPacketCompleteNat_intro
    (packet : ContactParameterCalibrationPacketNat)
    (hpairs : packet.contactPairRecords > 0)
    (hparameters : packet.parameterRecords >= packet.contactPairRecords)
    (hcolumns : packet.measurementColumns > 0)
    (hfit : packet.fitTemplateRows >= packet.contactPairRecords)
    (hholdout : packet.holdoutTemplateRows >= packet.contactPairRecords)
    (hmanifest : packet.artifactManifestEntries >= 5)
    (hprofile : packet.profileEvidenceReady > 0)
    (hmissing : packet.missingEvidence = 0) :
    contactParameterCalibrationPacketCompleteNat packet := by
  unfold contactParameterCalibrationPacketCompleteNat
  exact And.intro hpairs
    (And.intro hparameters
      (And.intro hcolumns
        (And.intro hfit
          (And.intro hholdout
            (And.intro hmanifest
              (And.intro hprofile hmissing))))))

theorem contactParameterCalibrationPacketCompleteNat_has_pairs
    {packet : ContactParameterCalibrationPacketNat}
    (hcomplete : contactParameterCalibrationPacketCompleteNat packet) :
    packet.contactPairRecords > 0 := by
  exact hcomplete.1

theorem contactParameterCalibrationPacketCompleteNat_parameters_cover_pairs
    {packet : ContactParameterCalibrationPacketNat}
    (hcomplete : contactParameterCalibrationPacketCompleteNat packet) :
    packet.parameterRecords >= packet.contactPairRecords := by
  exact hcomplete.2.1

theorem contactParameterCalibrationPacketCompleteNat_has_measurement_columns
    {packet : ContactParameterCalibrationPacketNat}
    (hcomplete : contactParameterCalibrationPacketCompleteNat packet) :
    packet.measurementColumns > 0 := by
  exact hcomplete.2.2.1

theorem contactParameterCalibrationPacketCompleteNat_fit_rows_cover_pairs
    {packet : ContactParameterCalibrationPacketNat}
    (hcomplete : contactParameterCalibrationPacketCompleteNat packet) :
    packet.fitTemplateRows >= packet.contactPairRecords := by
  exact hcomplete.2.2.2.1

theorem contactParameterCalibrationPacketCompleteNat_holdout_rows_cover_pairs
    {packet : ContactParameterCalibrationPacketNat}
    (hcomplete : contactParameterCalibrationPacketCompleteNat packet) :
    packet.holdoutTemplateRows >= packet.contactPairRecords := by
  exact hcomplete.2.2.2.2.1

theorem contactParameterCalibrationPacketCompleteNat_has_manifest
    {packet : ContactParameterCalibrationPacketNat}
    (hcomplete : contactParameterCalibrationPacketCompleteNat packet) :
    packet.artifactManifestEntries >= 5 := by
  exact hcomplete.2.2.2.2.2.1

theorem contactParameterCalibrationPacketCompleteNat_has_profile_evidence
    {packet : ContactParameterCalibrationPacketNat}
    (hcomplete : contactParameterCalibrationPacketCompleteNat packet) :
    packet.profileEvidenceReady > 0 := by
  exact hcomplete.2.2.2.2.2.2.1

theorem contactParameterCalibrationPacketCompleteNat_zero_missing_evidence
    {packet : ContactParameterCalibrationPacketNat}
    (hcomplete : contactParameterCalibrationPacketCompleteNat packet) :
    packet.missingEvidence = 0 := by
  exact hcomplete.2.2.2.2.2.2.2

theorem contactParameterBenchValidationReadyNat_intro
    (validation : ContactParameterBenchValidationNat)
    (hpacket : validation.packetReady > 0)
    (hfit : validation.fitRows > 0)
    (hholdout : validation.holdoutRows > 0)
    (hmeasurements : validation.completedMeasurements >= validation.fitRows + validation.holdoutRows)
    (hresiduals : validation.parameterResiduals > 0)
    (hfitpass : validation.fitPass > 0)
    (hholdoutpass : validation.holdoutPass > 0)
    (hindependent : validation.independentHoldoutPass > 0)
    (hmissing : validation.missingEvidence = 0) :
    contactParameterBenchValidationReadyNat validation := by
  unfold contactParameterBenchValidationReadyNat
  exact And.intro hpacket
    (And.intro hfit
      (And.intro hholdout
        (And.intro hmeasurements
          (And.intro hresiduals
            (And.intro hfitpass
              (And.intro hholdoutpass
                (And.intro hindependent hmissing)))))))

theorem contactParameterBenchValidationReadyNat_has_packet
    {validation : ContactParameterBenchValidationNat}
    (hready : contactParameterBenchValidationReadyNat validation) :
    validation.packetReady > 0 := by
  exact hready.1

theorem contactParameterBenchValidationReadyNat_has_fit_rows
    {validation : ContactParameterBenchValidationNat}
    (hready : contactParameterBenchValidationReadyNat validation) :
    validation.fitRows > 0 := by
  exact hready.2.1

theorem contactParameterBenchValidationReadyNat_has_holdout_rows
    {validation : ContactParameterBenchValidationNat}
    (hready : contactParameterBenchValidationReadyNat validation) :
    validation.holdoutRows > 0 := by
  exact hready.2.2.1

theorem contactParameterBenchValidationReadyNat_measurements_cover_rows
    {validation : ContactParameterBenchValidationNat}
    (hready : contactParameterBenchValidationReadyNat validation) :
    validation.completedMeasurements >= validation.fitRows + validation.holdoutRows := by
  exact hready.2.2.2.1

theorem contactParameterBenchValidationReadyNat_has_residuals
    {validation : ContactParameterBenchValidationNat}
    (hready : contactParameterBenchValidationReadyNat validation) :
    validation.parameterResiduals > 0 := by
  exact hready.2.2.2.2.1

theorem contactParameterBenchValidationReadyNat_has_fit_pass
    {validation : ContactParameterBenchValidationNat}
    (hready : contactParameterBenchValidationReadyNat validation) :
    validation.fitPass > 0 := by
  exact hready.2.2.2.2.2.1

theorem contactParameterBenchValidationReadyNat_has_holdout_pass
    {validation : ContactParameterBenchValidationNat}
    (hready : contactParameterBenchValidationReadyNat validation) :
    validation.holdoutPass > 0 := by
  exact hready.2.2.2.2.2.2.1

theorem contactParameterBenchValidationReadyNat_has_independent_holdout
    {validation : ContactParameterBenchValidationNat}
    (hready : contactParameterBenchValidationReadyNat validation) :
    validation.independentHoldoutPass > 0 := by
  exact hready.2.2.2.2.2.2.2.1

theorem contactParameterBenchValidationReadyNat_zero_missing_evidence
    {validation : ContactParameterBenchValidationNat}
    (hready : contactParameterBenchValidationReadyNat validation) :
    validation.missingEvidence = 0 := by
  exact hready.2.2.2.2.2.2.2.2

theorem contactParameterIntervalCalibrationReadyNat_intro
    (calibration : ContactParameterIntervalCalibrationNat)
    (hvalidation : calibration.benchValidationReady > 0)
    (hintervals : calibration.parameterIntervals > 0)
    (haccepted : calibration.acceptedParameterIntervals >= calibration.parameterIntervals)
    (huncertainty : calibration.uncertaintyRecords >= calibration.parameterIntervals)
    (hholdout : calibration.holdoutAgreementRecords >= calibration.parameterIntervals)
    (hinside : calibration.simulatorParametersInsideBounds > 0)
    (hmissing : calibration.missingEvidence = 0) :
    contactParameterIntervalCalibrationReadyNat calibration := by
  unfold contactParameterIntervalCalibrationReadyNat
  exact And.intro hvalidation
    (And.intro hintervals
      (And.intro haccepted
        (And.intro huncertainty
          (And.intro hholdout
            (And.intro hinside hmissing)))))

theorem contactParameterIntervalCalibrationReadyNat_has_validation
    {calibration : ContactParameterIntervalCalibrationNat}
    (hready : contactParameterIntervalCalibrationReadyNat calibration) :
    calibration.benchValidationReady > 0 := by
  exact hready.1

theorem contactParameterIntervalCalibrationReadyNat_has_intervals
    {calibration : ContactParameterIntervalCalibrationNat}
    (hready : contactParameterIntervalCalibrationReadyNat calibration) :
    calibration.parameterIntervals > 0 := by
  exact hready.2.1

theorem contactParameterIntervalCalibrationReadyNat_accepts_all_intervals
    {calibration : ContactParameterIntervalCalibrationNat}
    (hready : contactParameterIntervalCalibrationReadyNat calibration) :
    calibration.acceptedParameterIntervals >= calibration.parameterIntervals := by
  exact hready.2.2.1

theorem contactParameterIntervalCalibrationReadyNat_has_uncertainty
    {calibration : ContactParameterIntervalCalibrationNat}
    (hready : contactParameterIntervalCalibrationReadyNat calibration) :
    calibration.uncertaintyRecords >= calibration.parameterIntervals := by
  exact hready.2.2.2.1

theorem contactParameterIntervalCalibrationReadyNat_has_holdout_agreement
    {calibration : ContactParameterIntervalCalibrationNat}
    (hready : contactParameterIntervalCalibrationReadyNat calibration) :
    calibration.holdoutAgreementRecords >= calibration.parameterIntervals := by
  exact hready.2.2.2.2.1

theorem contactParameterIntervalCalibrationReadyNat_parameters_inside_bounds
    {calibration : ContactParameterIntervalCalibrationNat}
    (hready : contactParameterIntervalCalibrationReadyNat calibration) :
    calibration.simulatorParametersInsideBounds > 0 := by
  exact hready.2.2.2.2.2.1

theorem contactParameterIntervalCalibrationReadyNat_zero_missing_evidence
    {calibration : ContactParameterIntervalCalibrationNat}
    (hready : contactParameterIntervalCalibrationReadyNat calibration) :
    calibration.missingEvidence = 0 := by
  exact hready.2.2.2.2.2.2

theorem externalPhysicsRunReadyNat_intro
    (run : ExternalPhysicsRunNat)
    (hengine : run.engineAvailable > 0)
    (hexport : run.modelExportReady > 0)
    (hbodies : run.bodyRecords > 0)
    (hresults : run.resultBodyRecords >= run.bodyRecords)
    (hsteps : run.stepsRequested > 0)
    (hcompleted : run.stepsCompleted >= run.stepsRequested)
    (hmissing : run.missingEvidence = 0) :
    externalPhysicsRunReadyNat run := by
  unfold externalPhysicsRunReadyNat
  exact And.intro hengine
    (And.intro hexport
      (And.intro hbodies
        (And.intro hresults
          (And.intro hsteps
            (And.intro hcompleted hmissing)))))

theorem externalPhysicsRunReadyNat_has_engine
    {run : ExternalPhysicsRunNat}
    (hready : externalPhysicsRunReadyNat run) :
    run.engineAvailable > 0 := by
  exact hready.1

theorem externalPhysicsRunReadyNat_has_export
    {run : ExternalPhysicsRunNat}
    (hready : externalPhysicsRunReadyNat run) :
    run.modelExportReady > 0 := by
  exact hready.2.1

theorem externalPhysicsRunReadyNat_has_bodies
    {run : ExternalPhysicsRunNat}
    (hready : externalPhysicsRunReadyNat run) :
    run.bodyRecords > 0 := by
  exact hready.2.2.1

theorem externalPhysicsRunReadyNat_results_cover_bodies
    {run : ExternalPhysicsRunNat}
    (hready : externalPhysicsRunReadyNat run) :
    run.resultBodyRecords >= run.bodyRecords := by
  exact hready.2.2.2.1

theorem externalPhysicsRunReadyNat_has_steps
    {run : ExternalPhysicsRunNat}
    (hready : externalPhysicsRunReadyNat run) :
    run.stepsRequested > 0 := by
  exact hready.2.2.2.2.1

theorem externalPhysicsRunReadyNat_steps_completed
    {run : ExternalPhysicsRunNat}
    (hready : externalPhysicsRunReadyNat run) :
    run.stepsCompleted >= run.stepsRequested := by
  exact hready.2.2.2.2.2.1

theorem externalPhysicsRunReadyNat_zero_missing_evidence
    {run : ExternalPhysicsRunNat}
    (hready : externalPhysicsRunReadyNat run) :
    run.missingEvidence = 0 := by
  exact hready.2.2.2.2.2.2

theorem externalPhysicsComparisonReadyNat_intro
    (comparison : ExternalPhysicsComparisonNat)
    (hrun : comparison.runReady > 0)
    (hsimulator : comparison.simulatorRecords > 0)
    (hexternal : comparison.externalRecords >= comparison.simulatorRecords)
    (hmatched : comparison.matchedRecords >= comparison.simulatorRecords)
    (htolerance : comparison.toleranceRecords > 0)
    (hwithin : comparison.withinTolerance > 0)
    (hmissing : comparison.missingEvidence = 0) :
    externalPhysicsComparisonReadyNat comparison := by
  unfold externalPhysicsComparisonReadyNat
  exact And.intro hrun
    (And.intro hsimulator
      (And.intro hexternal
        (And.intro hmatched
          (And.intro htolerance
            (And.intro hwithin hmissing)))))

theorem externalPhysicsComparisonReadyNat_has_run
    {comparison : ExternalPhysicsComparisonNat}
    (hready : externalPhysicsComparisonReadyNat comparison) :
    comparison.runReady > 0 := by
  exact hready.1

theorem externalPhysicsComparisonReadyNat_has_simulator_records
    {comparison : ExternalPhysicsComparisonNat}
    (hready : externalPhysicsComparisonReadyNat comparison) :
    comparison.simulatorRecords > 0 := by
  exact hready.2.1

theorem externalPhysicsComparisonReadyNat_external_covers_simulator
    {comparison : ExternalPhysicsComparisonNat}
    (hready : externalPhysicsComparisonReadyNat comparison) :
    comparison.externalRecords >= comparison.simulatorRecords := by
  exact hready.2.2.1

theorem externalPhysicsComparisonReadyNat_matched_covers_simulator
    {comparison : ExternalPhysicsComparisonNat}
    (hready : externalPhysicsComparisonReadyNat comparison) :
    comparison.matchedRecords >= comparison.simulatorRecords := by
  exact hready.2.2.2.1

theorem externalPhysicsComparisonReadyNat_has_tolerance
    {comparison : ExternalPhysicsComparisonNat}
    (hready : externalPhysicsComparisonReadyNat comparison) :
    comparison.toleranceRecords > 0 := by
  exact hready.2.2.2.2.1

theorem externalPhysicsComparisonReadyNat_within_tolerance
    {comparison : ExternalPhysicsComparisonNat}
    (hready : externalPhysicsComparisonReadyNat comparison) :
    comparison.withinTolerance > 0 := by
  exact hready.2.2.2.2.2.1

theorem externalPhysicsComparisonReadyNat_zero_missing_evidence
    {comparison : ExternalPhysicsComparisonNat}
    (hready : externalPhysicsComparisonReadyNat comparison) :
    comparison.missingEvidence = 0 := by
  exact hready.2.2.2.2.2.2

theorem equilibriumRelationReadyNat_intro
    (relation : EquilibriumRelationNat)
    (hrealization : relation.realizationReady > 0)
    (hsolver : relation.solverSuccess > 0)
    (hcontact : relation.contactStateReady > 0)
    (henergy : relation.energyTerms > 0)
    (hnonnegative : relation.nonnegativeEnergyTerms >= relation.energyTerms)
    (hresidual : relation.balanceResidual <= relation.residualTolerance)
    (hmissing : relation.missingEvidence = 0) :
    equilibriumRelationReadyNat relation := by
  unfold equilibriumRelationReadyNat
  exact And.intro hrealization
    (And.intro hsolver
      (And.intro hcontact
        (And.intro henergy
          (And.intro hnonnegative
            (And.intro hresidual hmissing)))))

theorem equilibriumRelationReadyNat_has_realization
    {relation : EquilibriumRelationNat}
    (hready : equilibriumRelationReadyNat relation) :
    relation.realizationReady > 0 := by
  exact hready.1

theorem equilibriumRelationReadyNat_has_solver
    {relation : EquilibriumRelationNat}
    (hready : equilibriumRelationReadyNat relation) :
    relation.solverSuccess > 0 := by
  exact hready.2.1

theorem equilibriumRelationReadyNat_has_contact
    {relation : EquilibriumRelationNat}
    (hready : equilibriumRelationReadyNat relation) :
    relation.contactStateReady > 0 := by
  exact hready.2.2.1

theorem equilibriumRelationReadyNat_has_energy_terms
    {relation : EquilibriumRelationNat}
    (hready : equilibriumRelationReadyNat relation) :
    relation.energyTerms > 0 := by
  exact hready.2.2.2.1

theorem equilibriumRelationReadyNat_energy_terms_nonnegative
    {relation : EquilibriumRelationNat}
    (hready : equilibriumRelationReadyNat relation) :
    relation.nonnegativeEnergyTerms >= relation.energyTerms := by
  exact hready.2.2.2.2.1

theorem equilibriumRelationReadyNat_residual_within_tolerance
    {relation : EquilibriumRelationNat}
    (hready : equilibriumRelationReadyNat relation) :
    relation.balanceResidual <= relation.residualTolerance := by
  exact hready.2.2.2.2.2.1

theorem equilibriumRelationReadyNat_zero_missing_evidence
    {relation : EquilibriumRelationNat}
    (hready : equilibriumRelationReadyNat relation) :
    relation.missingEvidence = 0 := by
  exact hready.2.2.2.2.2.2

theorem reachableEquilibriumControllabilityReadyNat_intro
    (reachability : ReachableEquilibriumControllabilityNat)
    (hequilibrium : reachability.equilibriumReady > 0)
    (hbasis : reachability.actuatorBasis > 0)
    (hcolumns : reachability.responseColumns > 0)
    (htargets : reachability.targetCells > 0)
    (hresponses : reachability.reachableResponses > 0)
    (htopology : reachability.topologyEvidence > 0)
    (htargetPolicy :
      reachability.fullTargetRequired = 0 ∨ reachability.targetUnderactuated = 0)
    (hmissing : reachability.missingEvidence = 0) :
    reachableEquilibriumControllabilityReadyNat reachability := by
  unfold reachableEquilibriumControllabilityReadyNat
  exact And.intro hequilibrium
    (And.intro hbasis
      (And.intro hcolumns
        (And.intro htargets
          (And.intro hresponses
            (And.intro htopology
              (And.intro htargetPolicy hmissing))))))

theorem reachableEquilibriumControllabilityReadyNat_has_equilibrium
    {reachability : ReachableEquilibriumControllabilityNat}
    (hready : reachableEquilibriumControllabilityReadyNat reachability) :
    reachability.equilibriumReady > 0 := by
  exact hready.1

theorem reachableEquilibriumControllabilityReadyNat_has_actuator_basis
    {reachability : ReachableEquilibriumControllabilityNat}
    (hready : reachableEquilibriumControllabilityReadyNat reachability) :
    reachability.actuatorBasis > 0 := by
  exact hready.2.1

theorem reachableEquilibriumControllabilityReadyNat_has_response_columns
    {reachability : ReachableEquilibriumControllabilityNat}
    (hready : reachableEquilibriumControllabilityReadyNat reachability) :
    reachability.responseColumns > 0 := by
  exact hready.2.2.1

theorem reachableEquilibriumControllabilityReadyNat_has_target_cells
    {reachability : ReachableEquilibriumControllabilityNat}
    (hready : reachableEquilibriumControllabilityReadyNat reachability) :
    reachability.targetCells > 0 := by
  exact hready.2.2.2.1

theorem reachableEquilibriumControllabilityReadyNat_has_reachable_responses
    {reachability : ReachableEquilibriumControllabilityNat}
    (hready : reachableEquilibriumControllabilityReadyNat reachability) :
    reachability.reachableResponses > 0 := by
  exact hready.2.2.2.2.1

theorem reachableEquilibriumControllabilityReadyNat_has_topology
    {reachability : ReachableEquilibriumControllabilityNat}
    (hready : reachableEquilibriumControllabilityReadyNat reachability) :
    reachability.topologyEvidence > 0 := by
  exact hready.2.2.2.2.2.1

theorem reachableEquilibriumControllabilityReadyNat_target_policy
    {reachability : ReachableEquilibriumControllabilityNat}
    (hready : reachableEquilibriumControllabilityReadyNat reachability) :
    reachability.fullTargetRequired = 0 ∨ reachability.targetUnderactuated = 0 := by
  exact hready.2.2.2.2.2.2.1

theorem reachableEquilibriumControllabilityReadyNat_zero_missing_evidence
    {reachability : ReachableEquilibriumControllabilityNat}
    (hready : reachableEquilibriumControllabilityReadyNat reachability) :
    reachability.missingEvidence = 0 := by
  exact hready.2.2.2.2.2.2.2

theorem reachableEquilibriumBenchProtocolReadyNat_intro
    (protocol : ReachableEquilibriumBenchProtocolNat)
    (hreachability : protocol.reachabilityReady > 0)
    (hsteps : protocol.scenarioCount > 0)
    (htrials : protocol.actuatorColumnTrials > 0)
    (htargets : protocol.targetObservationCells > 0)
    (hmeasurements : protocol.measurementColumns > 0)
    (htopology : protocol.topologyControlPolicy > 0)
    (hpassfail : protocol.passFailCriteria > 0)
    (hmissing : protocol.missingEvidence = 0) :
    reachableEquilibriumBenchProtocolReadyNat protocol := by
  unfold reachableEquilibriumBenchProtocolReadyNat
  exact And.intro hreachability
    (And.intro hsteps
      (And.intro htrials
        (And.intro htargets
          (And.intro hmeasurements
            (And.intro htopology
              (And.intro hpassfail hmissing))))))

theorem reachableEquilibriumBenchProtocolReadyNat_has_reachability
    {protocol : ReachableEquilibriumBenchProtocolNat}
    (hready : reachableEquilibriumBenchProtocolReadyNat protocol) :
    protocol.reachabilityReady > 0 := by
  exact hready.1

theorem reachableEquilibriumBenchProtocolReadyNat_has_steps
    {protocol : ReachableEquilibriumBenchProtocolNat}
    (hready : reachableEquilibriumBenchProtocolReadyNat protocol) :
    protocol.scenarioCount > 0 := by
  exact hready.2.1

theorem reachableEquilibriumBenchProtocolReadyNat_has_actuator_trials
    {protocol : ReachableEquilibriumBenchProtocolNat}
    (hready : reachableEquilibriumBenchProtocolReadyNat protocol) :
    protocol.actuatorColumnTrials > 0 := by
  exact hready.2.2.1

theorem reachableEquilibriumBenchProtocolReadyNat_has_targets
    {protocol : ReachableEquilibriumBenchProtocolNat}
    (hready : reachableEquilibriumBenchProtocolReadyNat protocol) :
    protocol.targetObservationCells > 0 := by
  exact hready.2.2.2.1

theorem reachableEquilibriumBenchProtocolReadyNat_has_measurement_columns
    {protocol : ReachableEquilibriumBenchProtocolNat}
    (hready : reachableEquilibriumBenchProtocolReadyNat protocol) :
    protocol.measurementColumns > 0 := by
  exact hready.2.2.2.2.1

theorem reachableEquilibriumBenchProtocolReadyNat_has_topology_policy
    {protocol : ReachableEquilibriumBenchProtocolNat}
    (hready : reachableEquilibriumBenchProtocolReadyNat protocol) :
    protocol.topologyControlPolicy > 0 := by
  exact hready.2.2.2.2.2.1

theorem reachableEquilibriumBenchProtocolReadyNat_has_pass_fail
    {protocol : ReachableEquilibriumBenchProtocolNat}
    (hready : reachableEquilibriumBenchProtocolReadyNat protocol) :
    protocol.passFailCriteria > 0 := by
  exact hready.2.2.2.2.2.2.1

theorem reachableEquilibriumBenchProtocolReadyNat_zero_missing_evidence
    {protocol : ReachableEquilibriumBenchProtocolNat}
    (hready : reachableEquilibriumBenchProtocolReadyNat protocol) :
    protocol.missingEvidence = 0 := by
  exact hready.2.2.2.2.2.2.2

theorem reachableEquilibriumBenchValidationReadyNat_intro
    (validation : ReachableEquilibriumBenchValidationNat)
    (hprotocol : validation.protocolReady > 0)
    (hrows : validation.resultRows > 0)
    (htargets : validation.targetMeasurements > 0)
    (hcompleted : validation.completedMeasurements >= validation.targetMeasurements)
    (hreachability : validation.reachabilityChecks > 0)
    (htopology : validation.topologyLeakageChecks > 0)
    (hgroup : validation.groupSequenceChecks > 0)
    (hpass : validation.comparisonPass > 0)
    (hmissing : validation.missingEvidence = 0) :
    reachableEquilibriumBenchValidationReadyNat validation := by
  unfold reachableEquilibriumBenchValidationReadyNat
  exact And.intro hprotocol
    (And.intro hrows
      (And.intro htargets
        (And.intro hcompleted
          (And.intro hreachability
            (And.intro htopology
              (And.intro hgroup
                (And.intro hpass hmissing)))))))

theorem reachableEquilibriumBenchValidationReadyNat_has_protocol
    {validation : ReachableEquilibriumBenchValidationNat}
    (hready : reachableEquilibriumBenchValidationReadyNat validation) :
    validation.protocolReady > 0 := by
  exact hready.1

theorem reachableEquilibriumBenchValidationReadyNat_has_rows
    {validation : ReachableEquilibriumBenchValidationNat}
    (hready : reachableEquilibriumBenchValidationReadyNat validation) :
    validation.resultRows > 0 := by
  exact hready.2.1

theorem reachableEquilibriumBenchValidationReadyNat_has_targets
    {validation : ReachableEquilibriumBenchValidationNat}
    (hready : reachableEquilibriumBenchValidationReadyNat validation) :
    validation.targetMeasurements > 0 := by
  exact hready.2.2.1

theorem reachableEquilibriumBenchValidationReadyNat_measurements_complete
    {validation : ReachableEquilibriumBenchValidationNat}
    (hready : reachableEquilibriumBenchValidationReadyNat validation) :
    validation.completedMeasurements >= validation.targetMeasurements := by
  exact hready.2.2.2.1

theorem reachableEquilibriumBenchValidationReadyNat_has_reachability_checks
    {validation : ReachableEquilibriumBenchValidationNat}
    (hready : reachableEquilibriumBenchValidationReadyNat validation) :
    validation.reachabilityChecks > 0 := by
  exact hready.2.2.2.2.1

theorem reachableEquilibriumBenchValidationReadyNat_has_topology_checks
    {validation : ReachableEquilibriumBenchValidationNat}
    (hready : reachableEquilibriumBenchValidationReadyNat validation) :
    validation.topologyLeakageChecks > 0 := by
  exact hready.2.2.2.2.2.1

theorem reachableEquilibriumBenchValidationReadyNat_has_group_checks
    {validation : ReachableEquilibriumBenchValidationNat}
    (hready : reachableEquilibriumBenchValidationReadyNat validation) :
    validation.groupSequenceChecks > 0 := by
  exact hready.2.2.2.2.2.2.1

theorem reachableEquilibriumBenchValidationReadyNat_has_pass
    {validation : ReachableEquilibriumBenchValidationNat}
    (hready : reachableEquilibriumBenchValidationReadyNat validation) :
    validation.comparisonPass > 0 := by
  exact hready.2.2.2.2.2.2.2.1

theorem reachableEquilibriumBenchValidationReadyNat_zero_missing_evidence
    {validation : ReachableEquilibriumBenchValidationNat}
    (hready : reachableEquilibriumBenchValidationReadyNat validation) :
    validation.missingEvidence = 0 := by
  exact hready.2.2.2.2.2.2.2.2

theorem reachableEquilibriumAmplitudeCalibrationReadyNat_intro
    (calibration : ReachableEquilibriumAmplitudeCalibrationNat)
    (hvalidation : calibration.benchValidationReady > 0)
    (hrepeats : calibration.repeatedTrialGroups > 0)
    (hamplitudes : calibration.amplitudeEstimates > 0)
    (hresiduals : calibration.residualFieldCells > 0)
    (huncertainty : calibration.uncertaintyBands > 0)
    (htopology : calibration.topologyLeakageBands > 0)
    (hgroup : calibration.groupSequenceResiduals > 0)
    (hmissing : calibration.missingEvidence = 0) :
    reachableEquilibriumAmplitudeCalibrationReadyNat calibration := by
  unfold reachableEquilibriumAmplitudeCalibrationReadyNat
  exact And.intro hvalidation
    (And.intro hrepeats
      (And.intro hamplitudes
        (And.intro hresiduals
          (And.intro huncertainty
            (And.intro htopology
              (And.intro hgroup hmissing))))))

theorem reachableEquilibriumAmplitudeCalibrationReadyNat_has_validation
    {calibration : ReachableEquilibriumAmplitudeCalibrationNat}
    (hready : reachableEquilibriumAmplitudeCalibrationReadyNat calibration) :
    calibration.benchValidationReady > 0 := by
  exact hready.1

theorem reachableEquilibriumAmplitudeCalibrationReadyNat_has_repeated_trials
    {calibration : ReachableEquilibriumAmplitudeCalibrationNat}
    (hready : reachableEquilibriumAmplitudeCalibrationReadyNat calibration) :
    calibration.repeatedTrialGroups > 0 := by
  exact hready.2.1

theorem reachableEquilibriumAmplitudeCalibrationReadyNat_has_amplitudes
    {calibration : ReachableEquilibriumAmplitudeCalibrationNat}
    (hready : reachableEquilibriumAmplitudeCalibrationReadyNat calibration) :
    calibration.amplitudeEstimates > 0 := by
  exact hready.2.2.1

theorem reachableEquilibriumAmplitudeCalibrationReadyNat_has_residual_field
    {calibration : ReachableEquilibriumAmplitudeCalibrationNat}
    (hready : reachableEquilibriumAmplitudeCalibrationReadyNat calibration) :
    calibration.residualFieldCells > 0 := by
  exact hready.2.2.2.1

theorem reachableEquilibriumAmplitudeCalibrationReadyNat_has_uncertainty
    {calibration : ReachableEquilibriumAmplitudeCalibrationNat}
    (hready : reachableEquilibriumAmplitudeCalibrationReadyNat calibration) :
    calibration.uncertaintyBands > 0 := by
  exact hready.2.2.2.2.1

theorem reachableEquilibriumAmplitudeCalibrationReadyNat_has_topology_bands
    {calibration : ReachableEquilibriumAmplitudeCalibrationNat}
    (hready : reachableEquilibriumAmplitudeCalibrationReadyNat calibration) :
    calibration.topologyLeakageBands > 0 := by
  exact hready.2.2.2.2.2.1

theorem reachableEquilibriumAmplitudeCalibrationReadyNat_has_group_residuals
    {calibration : ReachableEquilibriumAmplitudeCalibrationNat}
    (hready : reachableEquilibriumAmplitudeCalibrationReadyNat calibration) :
    calibration.groupSequenceResiduals > 0 := by
  exact hready.2.2.2.2.2.2.1

theorem reachableEquilibriumAmplitudeCalibrationReadyNat_zero_missing_evidence
    {calibration : ReachableEquilibriumAmplitudeCalibrationNat}
    (hready : reachableEquilibriumAmplitudeCalibrationReadyNat calibration) :
    calibration.missingEvidence = 0 := by
  exact hready.2.2.2.2.2.2.2

theorem reachableEquilibriumEmpiricalProfileReadyNat_intro
    (profile : ReachableEquilibriumEmpiricalProfileNat)
    (hamplitude : profile.amplitudeCalibrationReady > 0)
    (hbounded : profile.boundedProposals > 0)
    (hsafe : profile.safeProposals > 0)
    (htolerances : profile.toleranceRecords > 0)
    (huncertainty : profile.uncertaintyRecords > 0)
    (hholdout : profile.holdoutHooks > 0)
    (hmissing : profile.missingEvidence = 0) :
    reachableEquilibriumEmpiricalProfileReadyNat profile := by
  unfold reachableEquilibriumEmpiricalProfileReadyNat
  exact And.intro hamplitude
    (And.intro hbounded
      (And.intro hsafe
        (And.intro htolerances
          (And.intro huncertainty
            (And.intro hholdout hmissing)))))

theorem reachableEquilibriumEmpiricalProfileReadyNat_has_amplitude
    {profile : ReachableEquilibriumEmpiricalProfileNat}
    (hready : reachableEquilibriumEmpiricalProfileReadyNat profile) :
    profile.amplitudeCalibrationReady > 0 := by
  exact hready.1

theorem reachableEquilibriumEmpiricalProfileReadyNat_has_bounded_proposals
    {profile : ReachableEquilibriumEmpiricalProfileNat}
    (hready : reachableEquilibriumEmpiricalProfileReadyNat profile) :
    profile.boundedProposals > 0 := by
  exact hready.2.1

theorem reachableEquilibriumEmpiricalProfileReadyNat_has_safe_proposals
    {profile : ReachableEquilibriumEmpiricalProfileNat}
    (hready : reachableEquilibriumEmpiricalProfileReadyNat profile) :
    profile.safeProposals > 0 := by
  exact hready.2.2.1

theorem reachableEquilibriumEmpiricalProfileReadyNat_has_tolerances
    {profile : ReachableEquilibriumEmpiricalProfileNat}
    (hready : reachableEquilibriumEmpiricalProfileReadyNat profile) :
    profile.toleranceRecords > 0 := by
  exact hready.2.2.2.1

theorem reachableEquilibriumEmpiricalProfileReadyNat_has_uncertainty
    {profile : ReachableEquilibriumEmpiricalProfileNat}
    (hready : reachableEquilibriumEmpiricalProfileReadyNat profile) :
    profile.uncertaintyRecords > 0 := by
  exact hready.2.2.2.2.1

theorem reachableEquilibriumEmpiricalProfileReadyNat_has_holdout_hooks
    {profile : ReachableEquilibriumEmpiricalProfileNat}
    (hready : reachableEquilibriumEmpiricalProfileReadyNat profile) :
    profile.holdoutHooks > 0 := by
  exact hready.2.2.2.2.2.1

theorem reachableEquilibriumEmpiricalProfileReadyNat_zero_missing_evidence
    {profile : ReachableEquilibriumEmpiricalProfileNat}
    (hready : reachableEquilibriumEmpiricalProfileReadyNat profile) :
    profile.missingEvidence = 0 := by
  exact hready.2.2.2.2.2.2

theorem reachableEquilibriumProfileInverseReadyNat_intro
    (inverse : ReachableEquilibriumProfileInverseNat)
    (hprofile : inverse.empiricalProfileReady > 0)
    (hsafe : inverse.safeProfileProposals > 0)
    (htargets : inverse.targetCells > 0)
    (hsolve : inverse.inverseSolveSuccess > 0)
    (hresiduals : inverse.residualRecords > 0)
    (hscores : inverse.scoreTerms > 0)
    (hreadonly : inverse.readOnlyProfileUse > 0)
    (hmissing : inverse.missingEvidence = 0) :
    reachableEquilibriumProfileInverseReadyNat inverse := by
  unfold reachableEquilibriumProfileInverseReadyNat
  exact And.intro hprofile
    (And.intro hsafe
      (And.intro htargets
          (And.intro hsolve
            (And.intro hresiduals
              (And.intro hscores
              (And.intro hreadonly hmissing))))))

theorem reachableEquilibriumProfileInverseReadyNat_has_profile
    {inverse : ReachableEquilibriumProfileInverseNat}
    (hready : reachableEquilibriumProfileInverseReadyNat inverse) :
    inverse.empiricalProfileReady > 0 := by
  exact hready.1

theorem reachableEquilibriumProfileInverseReadyNat_has_safe_proposals
    {inverse : ReachableEquilibriumProfileInverseNat}
    (hready : reachableEquilibriumProfileInverseReadyNat inverse) :
    inverse.safeProfileProposals > 0 := by
  exact hready.2.1

theorem reachableEquilibriumProfileInverseReadyNat_has_targets
    {inverse : ReachableEquilibriumProfileInverseNat}
    (hready : reachableEquilibriumProfileInverseReadyNat inverse) :
    inverse.targetCells > 0 := by
  exact hready.2.2.1

theorem reachableEquilibriumProfileInverseReadyNat_has_solve
    {inverse : ReachableEquilibriumProfileInverseNat}
    (hready : reachableEquilibriumProfileInverseReadyNat inverse) :
    inverse.inverseSolveSuccess > 0 := by
  exact hready.2.2.2.1

theorem reachableEquilibriumProfileInverseReadyNat_has_residuals
    {inverse : ReachableEquilibriumProfileInverseNat}
    (hready : reachableEquilibriumProfileInverseReadyNat inverse) :
    inverse.residualRecords > 0 := by
  exact hready.2.2.2.2.1

theorem reachableEquilibriumProfileInverseReadyNat_has_scores
    {inverse : ReachableEquilibriumProfileInverseNat}
    (hready : reachableEquilibriumProfileInverseReadyNat inverse) :
    inverse.scoreTerms > 0 := by
  exact hready.2.2.2.2.2.1

theorem reachableEquilibriumProfileInverseReadyNat_read_only_profile_use
    {inverse : ReachableEquilibriumProfileInverseNat}
    (hready : reachableEquilibriumProfileInverseReadyNat inverse) :
    inverse.readOnlyProfileUse > 0 := by
  exact hready.2.2.2.2.2.2.1

theorem reachableEquilibriumProfileInverseReadyNat_zero_missing_evidence
    {inverse : ReachableEquilibriumProfileInverseNat}
    (hready : reachableEquilibriumProfileInverseReadyNat inverse) :
    inverse.missingEvidence = 0 := by
  exact hready.2.2.2.2.2.2.2

theorem reachableEquilibriumProfileInverseAcceptanceReadyNat_intro
    (acceptance : ReachableEquilibriumProfileInverseAcceptanceNat)
    (hprofile : acceptance.profileInverseReady > 0)
    (hscoreRecords : acceptance.scoreRecords > 0)
    (hscore : acceptance.residualScore <= acceptance.residualLimit)
    (hbands : acceptance.bandFailures <= acceptance.maxBandFailures)
    (hactuators : acceptance.actuatorCount <= acceptance.actuatorLimit)
    (haccepted : acceptance.acceptedForPreview > 0)
    (hmissing : acceptance.missingEvidence = 0) :
    reachableEquilibriumProfileInverseAcceptanceReadyNat acceptance := by
  unfold reachableEquilibriumProfileInverseAcceptanceReadyNat
  exact And.intro hprofile
    (And.intro hscoreRecords
      (And.intro hscore
        (And.intro hbands
          (And.intro hactuators
            (And.intro haccepted hmissing)))))

theorem reachableEquilibriumProfileInverseAcceptanceReadyNat_has_profile_inverse
    {acceptance : ReachableEquilibriumProfileInverseAcceptanceNat}
    (hready : reachableEquilibriumProfileInverseAcceptanceReadyNat acceptance) :
    acceptance.profileInverseReady > 0 := by
  exact hready.1

theorem reachableEquilibriumProfileInverseAcceptanceReadyNat_has_score_records
    {acceptance : ReachableEquilibriumProfileInverseAcceptanceNat}
    (hready : reachableEquilibriumProfileInverseAcceptanceReadyNat acceptance) :
    acceptance.scoreRecords > 0 := by
  exact hready.2.1

theorem reachableEquilibriumProfileInverseAcceptanceReadyNat_score_within_limit
    {acceptance : ReachableEquilibriumProfileInverseAcceptanceNat}
    (hready : reachableEquilibriumProfileInverseAcceptanceReadyNat acceptance) :
    acceptance.residualScore <= acceptance.residualLimit := by
  exact hready.2.2.1

theorem reachableEquilibriumProfileInverseAcceptanceReadyNat_band_failures_within_limit
    {acceptance : ReachableEquilibriumProfileInverseAcceptanceNat}
    (hready : reachableEquilibriumProfileInverseAcceptanceReadyNat acceptance) :
    acceptance.bandFailures <= acceptance.maxBandFailures := by
  exact hready.2.2.2.1

theorem reachableEquilibriumProfileInverseAcceptanceReadyNat_actuators_within_limit
    {acceptance : ReachableEquilibriumProfileInverseAcceptanceNat}
    (hready : reachableEquilibriumProfileInverseAcceptanceReadyNat acceptance) :
    acceptance.actuatorCount <= acceptance.actuatorLimit := by
  exact hready.2.2.2.2.1

theorem reachableEquilibriumProfileInverseAcceptanceReadyNat_accepted
    {acceptance : ReachableEquilibriumProfileInverseAcceptanceNat}
    (hready : reachableEquilibriumProfileInverseAcceptanceReadyNat acceptance) :
    acceptance.acceptedForPreview > 0 := by
  exact hready.2.2.2.2.2.1

theorem reachableEquilibriumProfileInverseAcceptanceReadyNat_zero_missing_evidence
    {acceptance : ReachableEquilibriumProfileInverseAcceptanceNat}
    (hready : reachableEquilibriumProfileInverseAcceptanceReadyNat acceptance) :
    acceptance.missingEvidence = 0 := by
  exact hready.2.2.2.2.2.2

theorem reachableEquilibriumProfileInversePreviewPacketReadyNat_intro
    (packet : ReachableEquilibriumProfileInversePreviewPacketNat)
    (hacceptance : packet.acceptanceReady > 0)
    (hcommands : packet.commandRecords > 0)
    (hevents : packet.eventRecords >= packet.commandRecords)
    (htargets : packet.targetRecords > 0)
    (hresiduals : packet.residualRecords > 0)
    (hreadonly : packet.readOnlyPacket > 0)
    (hmissing : packet.missingEvidence = 0) :
    reachableEquilibriumProfileInversePreviewPacketReadyNat packet := by
  unfold reachableEquilibriumProfileInversePreviewPacketReadyNat
  exact And.intro hacceptance
    (And.intro hcommands
      (And.intro hevents
        (And.intro htargets
          (And.intro hresiduals
            (And.intro hreadonly hmissing)))))

theorem reachableEquilibriumProfileInversePreviewPacketReadyNat_has_acceptance
    {packet : ReachableEquilibriumProfileInversePreviewPacketNat}
    (hready : reachableEquilibriumProfileInversePreviewPacketReadyNat packet) :
    packet.acceptanceReady > 0 := by
  exact hready.1

theorem reachableEquilibriumProfileInversePreviewPacketReadyNat_has_commands
    {packet : ReachableEquilibriumProfileInversePreviewPacketNat}
    (hready : reachableEquilibriumProfileInversePreviewPacketReadyNat packet) :
    packet.commandRecords > 0 := by
  exact hready.2.1

theorem reachableEquilibriumProfileInversePreviewPacketReadyNat_events_cover_commands
    {packet : ReachableEquilibriumProfileInversePreviewPacketNat}
    (hready : reachableEquilibriumProfileInversePreviewPacketReadyNat packet) :
    packet.eventRecords >= packet.commandRecords := by
  exact hready.2.2.1

theorem reachableEquilibriumProfileInversePreviewPacketReadyNat_has_targets
    {packet : ReachableEquilibriumProfileInversePreviewPacketNat}
    (hready : reachableEquilibriumProfileInversePreviewPacketReadyNat packet) :
    packet.targetRecords > 0 := by
  exact hready.2.2.2.1

theorem reachableEquilibriumProfileInversePreviewPacketReadyNat_has_residuals
    {packet : ReachableEquilibriumProfileInversePreviewPacketNat}
    (hready : reachableEquilibriumProfileInversePreviewPacketReadyNat packet) :
    packet.residualRecords > 0 := by
  exact hready.2.2.2.2.1

theorem reachableEquilibriumProfileInversePreviewPacketReadyNat_read_only
    {packet : ReachableEquilibriumProfileInversePreviewPacketNat}
    (hready : reachableEquilibriumProfileInversePreviewPacketReadyNat packet) :
    packet.readOnlyPacket > 0 := by
  exact hready.2.2.2.2.2.1

theorem reachableEquilibriumProfileInversePreviewPacketReadyNat_zero_missing_evidence
    {packet : ReachableEquilibriumProfileInversePreviewPacketNat}
    (hready : reachableEquilibriumProfileInversePreviewPacketReadyNat packet) :
    packet.missingEvidence = 0 := by
  exact hready.2.2.2.2.2.2

theorem reachableEquilibriumProfileInversePreviewReplayReadyNat_intro
    (replay : ReachableEquilibriumProfileInversePreviewReplayNat)
    (hpacket : replay.packetReady > 0)
    (hcommands : replay.commandRecords > 0)
    (hreplayed : replay.replayedCommands = replay.commandRecords)
    (hevents : replay.eventCoverage > 0)
    (hsim : replay.simulationRecords > 0)
    (htargets : replay.targetResidualRecords > 0)
    (hagrees : replay.residualAgreement > 0)
    (hreadonly : replay.readOnlyReplay > 0)
    (hmissing : replay.missingEvidence = 0) :
    reachableEquilibriumProfileInversePreviewReplayReadyNat replay := by
  unfold reachableEquilibriumProfileInversePreviewReplayReadyNat
  exact And.intro hpacket
    (And.intro hcommands
      (And.intro hreplayed
        (And.intro hevents
          (And.intro hsim
            (And.intro htargets
              (And.intro hagrees
                (And.intro hreadonly hmissing)))))))

theorem reachableEquilibriumProfileInversePreviewReplayReadyNat_has_packet
    {replay : ReachableEquilibriumProfileInversePreviewReplayNat}
    (hready : reachableEquilibriumProfileInversePreviewReplayReadyNat replay) :
    replay.packetReady > 0 := by
  exact hready.1

theorem reachableEquilibriumProfileInversePreviewReplayReadyNat_has_commands
    {replay : ReachableEquilibriumProfileInversePreviewReplayNat}
    (hready : reachableEquilibriumProfileInversePreviewReplayReadyNat replay) :
    replay.commandRecords > 0 := by
  exact hready.2.1

theorem reachableEquilibriumProfileInversePreviewReplayReadyNat_replays_all_commands
    {replay : ReachableEquilibriumProfileInversePreviewReplayNat}
    (hready : reachableEquilibriumProfileInversePreviewReplayReadyNat replay) :
    replay.replayedCommands = replay.commandRecords := by
  exact hready.2.2.1

theorem reachableEquilibriumProfileInversePreviewReplayReadyNat_has_event_coverage
    {replay : ReachableEquilibriumProfileInversePreviewReplayNat}
    (hready : reachableEquilibriumProfileInversePreviewReplayReadyNat replay) :
    replay.eventCoverage > 0 := by
  exact hready.2.2.2.1

theorem reachableEquilibriumProfileInversePreviewReplayReadyNat_has_simulation
    {replay : ReachableEquilibriumProfileInversePreviewReplayNat}
    (hready : reachableEquilibriumProfileInversePreviewReplayReadyNat replay) :
    replay.simulationRecords > 0 := by
  exact hready.2.2.2.2.1

theorem reachableEquilibriumProfileInversePreviewReplayReadyNat_has_target_residuals
    {replay : ReachableEquilibriumProfileInversePreviewReplayNat}
    (hready : reachableEquilibriumProfileInversePreviewReplayReadyNat replay) :
    replay.targetResidualRecords > 0 := by
  exact hready.2.2.2.2.2.1

theorem reachableEquilibriumProfileInversePreviewReplayReadyNat_has_residual_agreement
    {replay : ReachableEquilibriumProfileInversePreviewReplayNat}
    (hready : reachableEquilibriumProfileInversePreviewReplayReadyNat replay) :
    replay.residualAgreement > 0 := by
  exact hready.2.2.2.2.2.2.1

theorem reachableEquilibriumProfileInversePreviewReplayReadyNat_read_only
    {replay : ReachableEquilibriumProfileInversePreviewReplayNat}
    (hready : reachableEquilibriumProfileInversePreviewReplayReadyNat replay) :
    replay.readOnlyReplay > 0 := by
  exact hready.2.2.2.2.2.2.2.1

theorem reachableEquilibriumProfileInversePreviewReplayReadyNat_zero_missing_evidence
    {replay : ReachableEquilibriumProfileInversePreviewReplayNat}
    (hready : reachableEquilibriumProfileInversePreviewReplayReadyNat replay) :
    replay.missingEvidence = 0 := by
  exact hready.2.2.2.2.2.2.2.2

theorem reachableEquilibriumProfileInversePreviewPhysicalReadyNat_intro
    (physical : ReachableEquilibriumProfileInversePreviewPhysicalNat)
    (hreplay : physical.replayReady > 0)
    (hsolver : physical.physicalSolverSuccess > 0)
    (hcommands : physical.commandRecords > 0)
    (htargets : physical.targetResidualRecords > 0)
    (hmodel : physical.modelComparisonRecords > 0)
    (henergy : physical.energyRecords > 0)
    (hreadonly : physical.readOnlyPhysicalPreview > 0)
    (hmissing : physical.missingEvidence = 0) :
    reachableEquilibriumProfileInversePreviewPhysicalReadyNat physical := by
  unfold reachableEquilibriumProfileInversePreviewPhysicalReadyNat
  exact And.intro hreplay
    (And.intro hsolver
      (And.intro hcommands
        (And.intro htargets
          (And.intro hmodel
            (And.intro henergy
              (And.intro hreadonly hmissing))))))

theorem reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_replay
    {physical : ReachableEquilibriumProfileInversePreviewPhysicalNat}
    (hready : reachableEquilibriumProfileInversePreviewPhysicalReadyNat physical) :
    physical.replayReady > 0 := by
  exact hready.1

theorem reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_solver
    {physical : ReachableEquilibriumProfileInversePreviewPhysicalNat}
    (hready : reachableEquilibriumProfileInversePreviewPhysicalReadyNat physical) :
    physical.physicalSolverSuccess > 0 := by
  exact hready.2.1

theorem reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_commands
    {physical : ReachableEquilibriumProfileInversePreviewPhysicalNat}
    (hready : reachableEquilibriumProfileInversePreviewPhysicalReadyNat physical) :
    physical.commandRecords > 0 := by
  exact hready.2.2.1

theorem reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_targets
    {physical : ReachableEquilibriumProfileInversePreviewPhysicalNat}
    (hready : reachableEquilibriumProfileInversePreviewPhysicalReadyNat physical) :
    physical.targetResidualRecords > 0 := by
  exact hready.2.2.2.1

theorem reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_model_comparison
    {physical : ReachableEquilibriumProfileInversePreviewPhysicalNat}
    (hready : reachableEquilibriumProfileInversePreviewPhysicalReadyNat physical) :
    physical.modelComparisonRecords > 0 := by
  exact hready.2.2.2.2.1

theorem reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_energy
    {physical : ReachableEquilibriumProfileInversePreviewPhysicalNat}
    (hready : reachableEquilibriumProfileInversePreviewPhysicalReadyNat physical) :
    physical.energyRecords > 0 := by
  exact hready.2.2.2.2.2.1

theorem reachableEquilibriumProfileInversePreviewPhysicalReadyNat_read_only
    {physical : ReachableEquilibriumProfileInversePreviewPhysicalNat}
    (hready : reachableEquilibriumProfileInversePreviewPhysicalReadyNat physical) :
    physical.readOnlyPhysicalPreview > 0 := by
  exact hready.2.2.2.2.2.2.1

theorem reachableEquilibriumProfileInversePreviewPhysicalReadyNat_zero_missing_evidence
    {physical : ReachableEquilibriumProfileInversePreviewPhysicalNat}
    (hready : reachableEquilibriumProfileInversePreviewPhysicalReadyNat physical) :
    physical.missingEvidence = 0 := by
  exact hready.2.2.2.2.2.2.2

theorem measuredEnergyResidualTripleNat_nonnegative
    (measuredSigned simulatedSigned
      measuredMagnitude simulatedMagnitude
      measuredContact simulatedContact : Nat) :
    0 <= measuredEnergyResidualTripleNat
      measuredSigned simulatedSigned
      measuredMagnitude simulatedMagnitude
      measuredContact simulatedContact := by
  exact Nat.zero_le _

theorem measuredEnergyResidualTripleNat_zero_when_equal
    (signed magnitude contact : Nat) :
    measuredEnergyResidualTripleNat
      signed signed
      magnitude magnitude
      contact contact = 0 := by
  unfold measuredEnergyResidualTripleNat
  simp [measuredWorkResidualNat_zero_when_equal]

theorem signedEnergyResidualTripleInt_nonnegative
    (measuredSigned simulatedSigned
      measuredMagnitude simulatedMagnitude
      measuredContact simulatedContact : Int) :
    0 <= signedEnergyResidualTripleInt
      measuredSigned simulatedSigned
      measuredMagnitude simulatedMagnitude
      measuredContact simulatedContact := by
  exact Nat.zero_le _

theorem signedEnergyResidualTripleInt_zero_when_equal
    (signed magnitude contact : Int) :
    signedEnergyResidualTripleInt
      signed signed
      magnitude magnitude
      contact contact = 0 := by
  unfold signedEnergyResidualTripleInt
  simp [signedWorkResidualInt_zero_when_equal]

theorem scaledSignedEnergyResidualTripleInt_numerator_nonnegative
    (measuredSigned simulatedSigned
      measuredMagnitude simulatedMagnitude
      measuredContact simulatedContact : Int)
    (denominator : Nat) :
    0 <=
      (scaledSignedEnergyResidualTripleInt
        measuredSigned simulatedSigned
        measuredMagnitude simulatedMagnitude
        measuredContact simulatedContact
        denominator).numerator := by
  exact Nat.zero_le _

theorem scaledSignedEnergyResidualTripleInt_zero_when_equal
    (signed magnitude contact : Int) (denominator : Nat) :
    (scaledSignedEnergyResidualTripleInt
      signed signed
      magnitude magnitude
      contact contact
      denominator).numerator = 0 := by
  unfold scaledSignedEnergyResidualTripleInt
  simp [signedEnergyResidualTripleInt_zero_when_equal]

theorem verticalLoadScenarioPassNat_zero_errors
    (tolerance : Nat) :
    verticalLoadScenarioPassNat 0 0 0 0 tolerance = true := by
  unfold verticalLoadScenarioPassNat
  simp

theorem verticalLoadScenarioPassNat_false_when_missing
    {missingMeasurements signedError magnitudeError contactError tolerance : Nat}
    (hmissing : 0 < missingMeasurements) :
    verticalLoadScenarioPassNat
      missingMeasurements signedError magnitudeError contactError tolerance = false := by
  unfold verticalLoadScenarioPassNat
  simp [Nat.ne_of_gt hmissing]

structure RigidBodyContactState (Body Contact : Type u) where
  activeContact : Contact -> Bool
  bodyOfContact : Contact -> Body

structure PinHoleClearanceModel (Scalar : Type u) where
  pinRadius : Scalar
  holeRadius : Scalar
  clearance : Scalar

structure RealValuedEnergyTargets (Scalar : Type u)
    [OfNat Scalar 0] [LE Scalar] where
  springEnergy : Scalar -> Scalar -> Scalar
  hingeEnergy : Scalar -> Scalar -> Scalar
  contactPenalty : Scalar -> Scalar -> Scalar
  spring_nonnegative :
    forall stiffness displacement,
      0 <= stiffness -> 0 <= springEnergy stiffness displacement
  hinge_nonnegative :
    forall stiffness angleError,
      0 <= stiffness -> 0 <= hingeEnergy stiffness angleError
  contact_nonnegative :
    forall stiffness penetration,
      0 <= stiffness -> 0 <= contactPenalty stiffness penetration

structure BenchValidationTargets (Scalar : Type u)
    [OfNat Scalar 0] [LE Scalar] where
  residual : Scalar -> Scalar -> Scalar
  residual_nonnegative :
    forall measured simulated, 0 <= residual measured simulated
  residual_zero_when_equal :
    forall value, residual value value = 0

end RadDiscontinuities
