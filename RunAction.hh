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

private:
    G4Accumulable<G4double> fSumTime = 0.;
    G4Accumulable<G4double> fSumTime2 = 0.;
    G4Accumulable<G4double> fSumRecoilEnergy = 0.;
    G4Accumulable<G4double> fSumInteractionDepth = 0.;
    G4Accumulable<G4double> fSumStoppingDepth = 0.;
    G4Accumulable<G4int> fStopCount = 0;
};

#endif
