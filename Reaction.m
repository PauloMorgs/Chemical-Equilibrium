classdef Reaction
    %Reaction
    % 
    
    properties
        Compounds            % Array of Compound objects
        Ni                   % Array with stoichiometric coefficients
        NCompounds           % Number of compouds in the reaction
        K_298                % Reaction equilibrium constant at 298 K
    end
    
    methods
        
        function obj = Reaction(Compounds, Ni)
            %Reaction construct an instance of this class
            
            if VerifyConsistency(obj, Compounds, Ni) ~= true
                error(strcat('Reaction not consistent. Verify stoichiometric coefficients.'));
            end
            
            obj.Compounds = Compounds;
            obj.Ni = Ni;
            obj.NCompounds = length(Compounds);
            obj.K_298 = CalculateK_298(obj, Compounds, Ni);
        end
        
        function K_T = CalculateK(obj, T, heat_capacity)
            %CalculateK
            % Analytical solution for the equilibrium constant based on 
            % heat capacity (Cp) constants A, B, C, D, E. dCoef is the 
            % array with delta A, delta B, ..., delta E. Same derivation 
            % from Koretsky p. 575.
            % User can choose between constant heat capacity or variable
            % heat capacity (function of temperature).
            %
            % Really sensitive to errors of the database -> exp(P+dP)
            
            % In case the user has set heat_capacity as constant
            % 
            if ~exist('heat_capacity', 'var')
                heat_capacity = 'variable';
            end
            
            
            Tref = 298.15;
            R = 8.314;
            
            % Calculation dH at T = 298.15 K
            % 
            dH_298 = CalculatedH_298(obj);
            
            % Verifying user options
            % 
            if strcmp(heat_capacity, 'constant') == 0
                
                % Calculation of ?A, ?B, ?C, ?D, ?E
                % e.g. ?A is defined by ?(A_i * v_i)
                % 
                
                dCoef = zeros(5,1);                                               % List with ?A, ?B, ?C, ?D, ?E
                for i = 1 : length(dCoef)                                         % Loop over the number of Cp coef.
                    for j = 1 : length(obj.Compounds)                             % Loop over the number of compounds
                        dCoef(i) = dCoef(i) + obj.Compounds(j).Cp(i) * obj.Ni(j); % Sum the each Cp coef. multiplied to stoichometric coef.
                    end
                end
                
                % Designing each coefficient to a variable.
                % (Human-interpretability reasons)
                dA = dCoef(1);
                dB = dCoef(2);
                dC = dCoef(3);
                dD = dCoef(4);
                dE = dCoef(5);
                
                % Calculation of K
                % 
                
                % The chemical equilibrium equation was divided in two
                % minor parts (a and b). (Koretsky p. 577)
                a = -dH_298/R + dA * Tref - (dB * Tref ^ 2) / 2 + (dC * Tref ^ 3) / 3 - dD / Tref + (dE * Tref ^ 4) / 4;
                b = dA * log(T / Tref) + (dB * (T - Tref)) / 2 + (dC * (T ^ 2 - Tref ^ 2)) / 6 + (dD * ((1 / T ^ 2)-(1 / Tref ^ 2))) / 2 + (dE * (T ^ 3 - Tref ^ 3)) / 12;
                
                % Full analytical equation of K considering dH_rxn as f(T)
                K_T = exp(a * (1 / T - 1 / Tref) + b) * obj.K_298;
                
            else
                
                % Calculation of K
                % 
                
                %Full analytical equation of K considering dH_rxn = constant
                K_T = exp((-dH_298 / R) * (1 / T - 1 / Tref)) * obj.K_298;
                
            end
            
        end
        
    end
    
    methods (Access = private)
        
        function boolean = VerifyConsistency(obj, Compounds, Ni)
            %VerifyConsistency
            % This method verifies the consistency between the atoms and 
            % the stoichiometric coefficients of the chemical equation
            
            % Matrix with the number of atoms of C, H, O, N for each 
            % compound in the respective rows
            % 
            
            % Declare shape of the matrix
            M = zeros(length(Compounds), 4);
            
            % Fit atoms in rows of the matrix
            for i = 1 : length(Compounds)
                M(i, :) = Compounds(i).Atoms';
            end
            
            % Verifying consistency by the stoichiometric coefficients
            % 
            
            % Using stoichiometric coefficients
            V = M.*Ni';
            
            % The sum of rows and columns must be zero in order to be a 
            % valid chemical equation
            if sum(sum(V)) == 0
                boolean = true;
            else
                boolean = false;
            end            
        end
        
        function K_298 =  CalculateK_298(obj, Compounds, Ni)
            %CalculateK_298
            % This method calculates the equilibrium contant at the 
            % reference temperature (Tref = 298.15 K).
            % Koretsky p. 572
            
            R = 8.314;
            T = 298.15;
            
            dg_rxn_298 = 0;
            for i = 1 : length(Compounds)
                dg_rxn_298 = dg_rxn_298 + (Compounds(i).GoT) * Ni(i);
            end
            
            K_298 = exp(-dg_rxn_298/(R * T));
        end
        
        function dH_298 =  CalculatedH_298(obj)
            %CalculateK_298
            % This method calculates the standard change in enthalpy of 
            % reaction at the reference temperature (Tref = 298.15 K).
            % Koretsky p. 572
            
            dH_298 = 0;
            for i = 1 : length(obj.Compounds)
                dH_298 = dH_298 + (obj.Compounds(i).HoT) * obj.Ni(i);
            end
            
        end
        
    end
end

















