#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"

class G4Run;

class RunAction : public G4UserRunAction {
public:
    RunAction();
    ~RunAction() override = default;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;

    void AddStoppingData(G4double time, G4double initialEnergy, G4double interactionDepth, G4double stoppingDepth);

    // Record a recoil that exits the target (thin-target case).
    // exitEnergy/exitBeta are signed: positive = forward exit, negative = backward exit.
    // Only forward exits enter the transmitted averages.
    void AddExitData(G4double initEnergy, G4double initBeta, G4double exitEnergy, G4double exitBeta, G4bool forward);

private:
    G4Accumulable<G4double> fSumTime = 0.;
    G4Accumulable<G4double> fSumTime2 = 0.;
    G4Accumulable<G4double> fSumRecoilEnergy = 0.;
    G4Accumulable<G4double> fSumInteractionDepth = 0.;
    G4Accumulable<G4double> fSumStoppingDepth = 0.;
    G4Accumulable<G4int> fStopCount = 0;

    // Transmitted-recoil statistics (forward exits only, except the counters)
    G4Accumulable<G4int> fTransmitCount = 0;
    G4Accumulable<G4int> fBackCount = 0;
    G4Accumulable<G4double> fSumTInitE = 0.;
    G4Accumulable<G4double> fSumTInitE2 = 0.;
    G4Accumulable<G4double> fSumTExitE = 0.;
    G4Accumulable<G4double> fSumTExitE2 = 0.;
    G4Accumulable<G4double> fSumTInitBeta = 0.;
    G4Accumulable<G4double> fSumTExitBeta = 0.;
};

#endif
