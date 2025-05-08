function output = blood_flow_units(input,unit)
% blood_flow_units prints a given value of volumetric blood flow rate in three
% commonly used units.
%
% inputs
%   input   double, numeric value to be converted
%   unit    string, either 'SI', 'CGS', or 'clinical'
%
% note: The CGS unit is g per cm^3 per s. The clinical unit is mL per 100 mL per
% min, which represents mL of blood per 100 mL of tissue per minute. The SI unit
% is s^-1.

    % convert input to SI units
    if strcmp(unit,'SI')

        temp = input;

    elseif strcmp(unit,'CGS')

        temp = input / 1.05;

    elseif strcmp(unit,'clinical')

        temp = input / 100 / 60;

    end

    % convert temp to each unit
    SI = temp;
    CGS = temp * 1.05;
    clinical = temp * 100 * 60;

    % store results
    output.SI = SI;
    output.CGS = CGS;
    output.clinical = clinical;

    % define units
    SI_unit = sprintf('[s^-1]');
    CGS_unit = sprintf('[g cm^-3 s^-1]');
    clinical_unit = sprintf('[mL 100 mL^-1 min^-1]');
    
    % display results
    if strcmp(unit,'SI')
        fprintf('\nThe input %f %s is:\n',input,SI_unit)
        fprintf('  %f %s\n',output.SI,SI_unit)
        fprintf('  %f %s\n',output.CGS,CGS_unit)
        fprintf('  %f %s\n',output.clinical,clinical_unit)
    elseif strcmp(unit,'CGS')
        fprintf('\nThe input %f %s is:\n',input,CGS_unit)
        fprintf('  %f %s\n',output.SI,SI_unit)
        fprintf('  %f %s\n',output.CGS,CGS_unit)
        fprintf('  %f %s\n',output.clinical,clinical_unit)
    elseif strcmp(unit,'clinical')
        fprintf('\nThe input %f %s is:\n',input,clinical_unit)
        fprintf('  %f %s\n',output.SI,SI_unit)
        fprintf('  %f %s\n',output.CGS,CGS_unit)
        fprintf('  %f %s\n',output.clinical,clinical_unit)
    end
end
