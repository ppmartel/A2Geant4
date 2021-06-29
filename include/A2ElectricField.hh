#ifndef A2ElectricField_h
#define A2ElectricField_h 1

#include "G4ElectricField.hh"
#include "G4UniformElectricField.hh"

class G4FieldManager;
class G4ChordFinder;
class G4EquationOfMotion;
class G4Mag_EqRhs;
class G4EqMagElectricField;
class G4MagIntegratorStepper;
class G4MagInt_Driver;


class A2ElectricField : public G4ElectricField
{
public:
	A2ElectricField();
	~A2ElectricField();

	void SetStepperType (G4int i) {fStepperType = i; CreateStepper(); }
	void SetMinStep(G4double s) { fMinStep = s; }

	void SetFieldValue(G4double fieldValue);
	void UpdateIntegrator();
	G4ElectricField* Construct(G4double fieldStrength);
	void GetFieldValue(const G4double point[4], G4double *field) const;

protected:
	G4FieldManager* GetGlobalFieldManager();
	void CreateStepper();

private:
	G4double fMinStep;
	G4FieldManager* fFieldManager;
	G4ChordFinder* fChordFinder;
	G4EqMagElectricField* fEquation;
	G4ElectricField* fField;
	G4MagIntegratorStepper* fStepper;
	G4MagInt_Driver* fIntegrationDriver;
	G4int fStepperType;
};
#endif
