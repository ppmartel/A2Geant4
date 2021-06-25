#include "A2ElectricField.hh"

#include "G4UniformElectricField.hh"
#include "G4UniformMagField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

A2ElectricField::A2ElectricField(){
	//constructor: let's do this the A2 way and here set everything to NULL
	fField=NULL;
	fEquation=NULL;
	fFieldManager=NULL;
	fStepper=NULL;
	fIntegrationDriver=NULL;
	fChordFinder=NULL;
	fStepperType=4; //choose a stepper
	fMinStep=0.010*mm; //initialize minimum step length
//	return fField;
}

A2ElectricField::~A2ElectricField(){
	//destructor: delete some shit
	delete fField;
	delete fEquation;
	delete fFieldManager;
	delete fStepper;
	delete fIntegrationDriver;
	delete fChordFinder;
}

G4ElectricField* A2ElectricField::Construct(){
	//do things
	fField= new G4UniformElectricField(G4ThreeVector(0.0,0.0,2.0*kilovolt/cm));
	fEquation= new G4EqMagElectricField(fField);
	fFieldManager = GetGlobalFieldManager();
	UpdateIntegrator();
	return fField;
}

void A2ElectricField::UpdateIntegrator(){
	CreateStepper();
	fIntegrationDriver = new G4MagInt_Driver(fMinStep, fStepper, fStepper->GetNumberOfVariables());
	fChordFinder = new G4ChordFinder(fIntegrationDriver);
	fFieldManager->SetChordFinder(fChordFinder);
	fFieldManager->SetDetectorField(fField);
}

void A2ElectricField::CreateStepper(){
	const G4int nvar = 8;

  	auto oldStepper= fStepper;

  	switch ( fStepperType )
  	{
  	  case 0:
  	    fStepper = new G4ExplicitEuler( fEquation, nvar );
  	    G4cout<<"G4ExplicitEuler is calledS"<<G4endl;
  	    break;
  	  case 1:
  	    fStepper = new G4ImplicitEuler( fEquation, nvar );
  	    G4cout<<"G4ImplicitEuler is called"<<G4endl;
  	    break;
  	  case 2:
  	    fStepper = new G4SimpleRunge( fEquation, nvar );
  	    G4cout<<"G4SimpleRunge is called"<<G4endl;
  	    break;
  	  case 3:
  	    fStepper = new G4SimpleHeum( fEquation, nvar );
  	    G4cout<<"G4SimpleHeum is called"<<G4endl;
  	    break;
  	  case 4:
  	    fStepper = new G4ClassicalRK4( fEquation, nvar );
  	    G4cout<<"G4ClassicalRK4 is called"<<G4endl;
  	    break;
  	  case 5:
  	    fStepper = new G4CashKarpRKF45( fEquation, nvar );
  	    G4cout<<"G4CashKarpRKF45 is called"<<G4endl;
  	    break;
  	  case 7:
  	    fStepper = 0; // new G4HelixExplicitEuler( fEquation );
  	    G4cout<<"G4HelixExplicitEuler is not valid for Electric Field"<<G4endl;
  	    break;
  	  case 8:
  	    fStepper = 0; // new G4HelixImplicitEuler( fEquation );
  	    G4cout<<"G4HelixImplicitEuler is not valid for Electric Field"<<G4endl;
  	    break;
  	  case 9:
  	    fStepper = 0; // new G4HelixSimpleRunge( fEquation );
  	    G4cout<<"G4HelixSimpleRunge is not valid for Electric Field"<<G4endl;
  	    break;
  	  default:  /* fStepper = 0; // Older code */
  	    fStepper = new G4ClassicalRK4( fEquation, nvar );
  	    G4cout<<"G4ClassicalRK4 (default) is called"<<G4endl;
  	    break;
  	}

  	delete oldStepper;
  
  	fIntegrationDriver->RenewStepperAndAdjust(fStepper);
}

//manually set value in Z
void A2ElectricField::SetFieldValue(G4double fieldValue){
	G4ThreeVector fieldVector(0.0,0.0,fieldValue*kilovolt/cm);
	G4FieldManager *fieldManager = GetGlobalFieldManager();
	if (fieldVector != G4ThreeVector(0.,0.,0.)){
		fField = new G4UniformElectricField(fieldVector);
	} else {
		fField=0;
	}
	fieldManager->SetDetectorField(fField);
	fEquation->SetFieldObj(fField);
}

//needs to be here for class to work
void A2ElectricField::GetFieldValue(const G4double point[4], G4double *field) const{
	//reimplement from G4ElectricField
	fField->GetFieldValue(point, field);
}



G4FieldManager* A2ElectricField::GetGlobalFieldManager(){
	return G4TransportationManager::GetTransportationManager()->GetFieldManager();
}
