

#include "A2SteppingAction.hh"

#include "A2DetectorConstruction.hh"
#include "A2EventAction.hh"

#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4Track.hh"
#include "G4Gamma.hh"
#include "G4Proton.hh"
#include "G4OpticalPhoton.hh"
#include "G4SDManager.hh"
#include "G4SteppingManager.hh"
#include "G4RunManager.hh"
#include "CLHEP/Units/SystemOfUnits.h"

using namespace CLHEP;

////#include "G4RunManager.hh"



A2SteppingAction::A2SteppingAction(A2DetectorConstruction* det,
                                         A2EventAction* evt)
{
    detector = det;
    eventaction = evt;
}



A2SteppingAction::~A2SteppingAction()
{ }



void A2SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4Track* track = aStep->GetTrack();

  G4StepPoint* startPoint = aStep->GetPreStepPoint();
  G4StepPoint* endPoint   = aStep->GetPostStepPoint();

  G4String particleName = track->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();

  //G4cout << particleName << "\t" << pds->GetProcessName() << "\t" << startPoint->GetPhysicalVolume()->GetName()<< " to " << endPoint->GetPhysicalVolume()->GetName();
  //G4cout << "\t" << startPoint->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial()->GetMaterialPropertiesTable()->GetConstProperty("SCINTILLATIONYIELD");
  //G4cout << G4endl;
  
  //Need to add some stuff here for sure! Follow AHeT code
  if (particleName == "e-") {
	  //propogate through field
	  //and see if you hit the anode
//G4https://root.cern/manual/histograms/histo-trial.png 
  }

  if (particleName == "opticalphoton") {
      G4OpBoundaryProcess* fOpProcess;

      // Retrieve the status of the photon
      G4OpBoundaryProcessStatus theStatus = Undefined;

      G4ProcessManager* OpManager =
              G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

      if (OpManager) {
          G4int MAXofPostStepLoops =
                  OpManager->GetPostStepProcessVector()->entries();
          G4ProcessVector* fPostStepDoItVector =
                  OpManager->GetPostStepProcessVector(typeDoIt);

          for ( G4int i=0; i<MAXofPostStepLoops; i++) {
              G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
              fOpProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
              if (fOpProcess) { theStatus = fOpProcess->GetStatus(); break;}
          }
      }

      // Detected by a detector
      if (theStatus == Detection) {
          //G4cout << "Detected in " << startPoint->GetPhysicalVolume()->GetName() << "\t" << aStep->GetTotalEnergyDeposit()/eV << G4endl;

          // Check if the photon hits the detector and process the hit if it does
          if ( startPoint->GetPhysicalVolume()->GetLogicalVolume()->GetName() == "LogicSiPMT" ) {

              G4SDManager* SDman = G4SDManager::GetSDMpointer();
              A2SD* AHe3SD = (A2SD*)SDman->FindSensitiveDetector("AHe3SD");

              //if (AHe3SD) AHe3SD->ProcessHits_AHe3(aStep, NULL);

              // Stop Tracking when it hits the detector's surface
              //ResetCounters();
              track->SetTrackStatus(fStopAndKill);
          }
      }
  }

//   G4VPhysicalVolume* volume = track->GetVolume();
  
//   // collect energy and track length step by step
   G4double edep = aStep->GetTotalEnergyDeposit();
  
//   G4double stepl = 0.;
//track electrons long enough to let them propogate to anode
if(track->GetDefinition()->GetParticleName()==G4String("e-")){
	//G4cout<<"Tracking electron "<<aStep->GetPreStepPoint()->GetGlobalTime()/ns<<" "<<track->GetKineticEnergy()/MeV<<" "<< fpSteppingManager->GetfCurrentVolume()->GetName()<<G4endl;
	if(aStep->GetPreStepPoint()->GetGlobalTime()>1000*ms)track->SetTrackStatus(fStopAndKill);
}
//stop tracking after the trigger time
else if(aStep->GetPreStepPoint()->GetGlobalTime()>2*ms)track->SetTrackStatus(fStopAndKill);
//   if(track->GetDefinition()->GetParticleName()==G4String("pi0"))
//     {G4cout<<"Got a pi0 "<<aStep->GetPreStepPoint()->GetGlobalTime()/ns<<" "<<track->GetKineticEnergy()/MeV<<" "<< fpSteppingManager->GetfCurrentVolume()->GetName()<<G4endl;track->SetTrackStatus(fStopAndKill);}
//  if(track->GetDefinition()->GetParticleName()==G4String("pi+"))
//    {G4cout<<"Got a pi+ "<<aStep->GetPreStepPoint()->GetGlobalTime()/ns<<" "<<track->GetKineticEnergy()/MeV<<" "<< fpSteppingManager->GetfCurrentVolume()->GetName()<<G4endl;track->SetTrackStatus(fStopAndKill);}
//  if(track->GetDefinition()->GetParticleName()==G4String("mu+"))
//    {G4cout<<"Got a mu+ "<<aStep->GetPreStepPoint()->GetGlobalTime()/ns<<" "<<track->GetKineticEnergy()/MeV<<" "<< fpSteppingManager->GetfCurrentVolume()->GetName()<<G4endl;track->SetTrackStatus(fStopAndKill);}
 //  if(track->GetDefinition()->GetParticleName()==G4String("pi0"))
//     {track->SetTrackStatus(fStopAndKill);}
//  if(track->GetDefinition()->GetParticleName()==G4String("pi+"))
//    {track->SetTrackStatus(fStopAndKill);}
//  if(track->GetDefinition()->GetParticleName()==G4String("pi-"))
//    {track->SetTrackStatus(fStopAndKill);}
//  if(track->GetDefinition()->GetParticleName()==G4String("mu+"))
//    {track->SetTrackStatus(fStopAndKill);}
  //if (track->GetDefinition()->GetPDGCharge() != 0.)
  // stepl = aStep->GetStepLength();
  //if(track->GetDefinition()->GetParticleName()==G4Gamma::Gamma()->GetParticleName()){

  //G4cout<<track->GetDefinition()->GetParticleName()<< track->GetTrackID()<<" process " <<fpSteppingManager->GetfCurrentProcess()->GetProcessName()<<" in "<< fpSteppingManager->GetfCurrentVolume()->GetName()<<G4endl;
  // }
  //  if(fpSteppingManager->GetfCurrentVolume()->GetName()==G4String("ANOIP"))G4cout<<"OK "<<track->GetDefinition()->GetParticleName()<<G4endl;
  //if(track->GetDefinition()->GetParticleName()==G4Proton::Proton()->GetParticleName()&&!(fpSteppingManager->GetfCurrentProcess()->GetProcessName()==G4String("msc"))&&!(fpSteppingManager->GetfCurrentProcess()->GetProcessName()==G4String("hIoni"))){
  // if(fpSteppingManager->GetfCurrentProcess()->GetProcessName()==G4String("msc"))
  // if(track->GetTrackID()==1&&fpSteppingManager->GetfCurrentVolume()->GetName()!=G4String("World")) G4cout<<track->GetDefinition()->GetParticleName()<< track->GetTrackID()<< " "<<track->GetParentID()<< " "<<track->GetKineticEnergy()<<" process " <<fpSteppingManager->GetfCurrentProcess()->GetProcessName()<<" in "<< fpSteppingManager->GetfCurrentVolume()->GetName()<<G4endl;
    //G4cout<<"Secondaries "<<aStep->GetSecondary()->size()<<" "<<aStep->GetfSecondary()<<G4endl;
  //}
  //G4cout <<" STEPPING ACTION "<<eventaction->GetNEvent()
  //  if(eventaction->GetNEvent()==1317){
  // G4cout<<track->GetDefinition()->GetParticleName()<< track->GetTrackID()<<" process " <<fpSteppingManager->GetfCurrentProcess()->GetProcessName()<<" in "<< fpSteppingManager->GetfCurrentVolume()->GetName()<< " "<<track->GetKineticEnergy()/MeV<<G4endl;
  // }

  //bug in phot process, can't get rid of gamma with energy 1.2E-5MeV
  //goes into infinite loop!
  if(track->GetDefinition()->GetParticleName()==G4Gamma::Gamma()->GetParticleName()&&track->GetKineticEnergy()/MeV<1E-4&&fpSteppingManager->GetfCurrentProcess()->GetProcessName()==G4String("phot"))track->SetTrackStatus(fStopAndKill);
}





