# InfantDENV

Code and data associated with the analysis in "Maternally-derived antibody titer dynamics and risk of hospitalised infant dengue disease", https://doi.org/10.1101/2022.11.18.22282500

### **Models of maternal antibody decay dynamics** 
```FitAbDecayModels.R``` calls the data and Stan models to reconstruct the decay dynamics of maternally-derived DENV antibodies.

### **Modelling titer-related mechanisms of infant disease risk**
```FitTiterRiskMechanisms.R``` calls data and Stan models for fitting titer-related mechanisms of infant disease risk to age-specific hospitalised infant dengue case data.  

### **Simulating long-term trends in age-specific infant disease risk**
```SimulateTemporalRisks.R``` calls the data and functions to simulate age-specific infant dengue disease risk trends over a 40 year period, for both the ADE (lognormal functional form) and non-ADE threshold scenarios. 
