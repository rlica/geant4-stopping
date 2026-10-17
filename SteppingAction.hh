#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4Step;

class SteppingAction : public G4UserSteppingAction {
public:
    SteppingAction() = default;
    ~SteppingAction() override = default;

    void UserSteppingAction(const G4Step* step) override;

private:
    // "Before the medium" state of the current recoil track, captured at its first step.
    G4double fRecoilInitEnergy = 0.0;
    G4double fRecoilInitBeta = 0.0;
};

#endif
