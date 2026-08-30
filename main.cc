#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "SteppingAction.hh"
#include "RunAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4EmParameters.hh"
#include "G4HadronicProcessStore.hh"

// Run with: razvan@Legion:~/geant4-tests/geant4-stopping/build$ make clean && make && ./stopping ../run.mac

// Custom minimal physics list for low-energy EM interactions only
class MinimalPhysicsList : public G4VModularPhysicsList {
public:
    MinimalPhysicsList() {
        SetVerboseLevel(0);
        RegisterPhysics(new G4EmStandardPhysics_option4(0)); 
        
        // Include inelastic and elastic nuclear interactions
        RegisterPhysics(new G4IonPhysics(0));
        RegisterPhysics(new G4HadronElasticPhysics(0));
        
        // Include Step Limiter for finer integration steps
        RegisterPhysics(new G4StepLimiterPhysics());

        // Tune EM parameters for high-accuracy low-energy tracking
        G4EmParameters* emParams = G4EmParameters::Instance();
        emParams->SetStepFunction(0.1, 0.01 * um); // tighter step constraints
        emParams->SetFluo(true); // Enable atomic deexcitation (fluorescence)
        emParams->SetAuger(true); // Enable Auger electrons
    }
};

// Minimal ActionInitialization
class ActionInitialization : public G4VUserActionInitialization {
public:
    ActionInitialization() = default;
    ~ActionInitialization() override = default;
    void BuildForMaster() const override {
        SetUserAction(new RunAction());
    }
    void Build() const override {
        SetUserAction(new PrimaryGeneratorAction());
        SetUserAction(new RunAction());
        SetUserAction(new SteppingAction());
    }
};

int main(int argc, char** argv) {
    // Force Serial RunManager for clean output
    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
    
    if (!runManager) {
        G4cerr << "ERR_MAIN_01: Failed to create G4RunManager." << G4endl;
        return 1;
    }

    // 1. Detector
    runManager->SetUserInitialization(new DetectorConstruction());

    // 2. Minimal Physics List
    G4VModularPhysicsList* physicsList = new MinimalPhysicsList();
    if (!physicsList) {
        G4cerr << "ERR_MAIN_02: Failed to allocate MinimalPhysicsList." << G4endl;
        return 1;
    }
    runManager->SetUserInitialization(physicsList);

    // 3. User Actions
    runManager->SetUserInitialization(new ActionInitialization());
    
    // Get UI pointer and silence all standard outputs before initialization
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    if (!UImanager) {
        G4cerr << "ERR_MAIN_03: Failed to get G4UImanager." << G4endl;
        return 1;
    }

    UImanager->ApplyCommand("/control/verbose 0");
    UImanager->ApplyCommand("/run/verbose 0");
    UImanager->ApplyCommand("/event/verbose 0");
    UImanager->ApplyCommand("/tracking/verbose 0");
    UImanager->ApplyCommand("/process/em/verbose 0");
    UImanager->ApplyCommand("/process/had/verbose 0");
    UImanager->ApplyCommand("/run/particle/verbose 0");

    // Force strict verbosity to 0 before initialization
    G4HadronicProcessStore::Instance()->SetVerbose(0);
    physicsList->SetVerboseLevel(0);

    // Initialize G4 kernel (moved to run.mac)
    // runManager->Initialize();

    // Execution
    if (argc == 1) {
        G4cout << "No macro provided. Running 10 test events..." << G4endl;
        UImanager->ApplyCommand("/run/beamOn 10");
    } else {
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    }

    delete runManager;
    return 0;
}
