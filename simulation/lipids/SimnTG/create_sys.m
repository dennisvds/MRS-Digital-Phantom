function sys = create_sys(FAchainlen, DBpos)
% The fatty acid chain length can be specified and 
% Chemical shifts and J coupling data were extracted from the following publications:
% For shifts and J coupling data see:
% K. Schaumburg and H._ J. Bernstein Lipids Vol3 1968 p 193
% Eleni Alexandri et al Molecules 2017, 22, 1663
% Carmen Salinero et al. Molecules 2012, 17, 6716-6727
% http://neu.ilps.org/wp-content/uploads/2019/07/Omega_3.pdf
%
% CH2 close to carboxyl are shifted
% list of chemical shift increments from: Frost and Gunstone, Chemistry and Physics of Lipids 15 (1975) 53-85
% position  delta   shift
% alpha     0.940	2.195
% beta      0.320	1.575
% gamma     0.060	1.315
% delta     0.000	1.255
% epsilon	0.025	1.280
% zeta      0.000	1.255
%
% CH2 shifts from double bond
% CH2 close to HC=CH		
% alpha     0.730	1.985
% beta      0.065	1.320
% gamma     0.020	1.275
% delta     0.015	1.270
% epsilon	0.005	1.260
% zeta      0.000	1.255
%
%
%%

%  Ronald Ouwerkerk NIDDK/NIH 2020


%% Set up spectrum paramaters

if nargin < 1 
    FAchainlen = 18;
end

% default to linoleic acid
if nargin < 2 
    DBpos = [6,9];
end

if 1 %~exist('TGspinSystems.mat') == 2
    %% Construct the resonances of the three FA and TG 
    sys = []; 
    %% The chain CH2 
                      
    % Set the chemical shift deltas due to shielding from the carboxyl (esther)
    carboxylshifts = [  0.940  % alpha
                        0.320  % beta
                        0.060  % gamma
                        0.000  % delta 
                        0.025  % epsilon
                        0.000];% zeta 

    % Set up shifts on CH2 caused by a double bond 
    unsatshifts = [ 0.730   % alpha
                    0.065   % beta
                    0.020   % gamma
                    0.015   % delta 
                    0.005]; % epsilon
                
    fppmCH3 =0.885;% base shift of the methyl group 0.885 ppm        
            
    % the shifts of the CH3 group due to double bond proximity          
    omegashiftsCH3 = [ fppmCH3  % Omega 7 and beyond
                       0.890    % Omega 6
                       0.905    % Omega 5
                       0.910    % Omega 4
                       0.955 ]; % Omega 3
                 
    %% Apply the chemical shift increments according to position in the FA chain
    % FA chain length includes COOH and CH3

    CH2s = 1:(FAchainlen-2); % create an array of positions
    fppmCH2 = 1.255; % The base shift for a FA chain methylene in ppm

    for ii = 1:6
        chainFCH2( ii ) = fppmCH2 + carboxylshifts(ii);
    end

    % the rest has the base shift
    for ii = 7:(FAchainlen-2)
        chainFCH2( ii ) = fppmCH2;
    end

    % The last one before the CH3 has a small shift and it is coupled with CH3
    chainFCH2( FAchainlen-2 ) = chainFCH2( FAchainlen-2 )+0.03;

    %%% Find double bond positions and remove them from the FA cahin CH2 list

    % Now lets count the other way: from CH3 down to COOH
    % Apply the chemical shift increments according to position in the FA chain
    % realtive to double bond
    % Identify the positions of the double bond carbons
    CH3toCOOpos = (flip( CH2s)+1);

    DBcarbons = false( 1, FAchainlen-2);
    alphaDBcarbons = zeros( 1, FAchainlen-2);
    
    sys = [];
    
    if ~isempty( DBpos  )
        for ii = 1: length( DBpos )
            DBidx = find( CH3toCOOpos == DBpos(ii));   
            % The other double bond is one closer to the carboxyl
            DBcarbons( [DBidx, DBidx-1] ) = true;
            % the same for adjacent CH2, the diallylic get a value 2
            alphaDBcarbons( [DBidx+1, DBidx-2] ) = alphaDBcarbons( [DBidx+1, DBidx-2] )+1;
        end

        NDB = length( find( DBcarbons ))/2;

        % Go through the FA chain carbons and adjust chemical shift for distance
        % from the double bond
        % Start from the COO- side
        do_CH2B = (DBidx-length(unsatshifts)):(DBidx-1);
        % fprintf(1, 'Adjusting shift for distance from double bond 1: from COO side\n');
        for ii = do_CH2B
            distance = -( ii -DBidx );
            % fprintf(1, '#%d distance=%d, shift %5.3f \n', ii, distance, unsatshifts(distance));
            chainFCH2( ii ) = chainFCH2( ii )+ unsatshifts(distance);
        end

        % Now through the FA chain carbons and adjust chemical shift for distance
        % from the double bond start from the double bond toward the CH3
        DBidx = find( DBcarbons, 1, 'last');
        do_CH2B = (DBidx+1):min( (FAchainlen-2), (DBidx+length(unsatshifts)));
        % fprintf(1, 'Adjusting shift for distance from double bond 2: from CH3 side\n');
        for ii = do_CH2B
            distance = ( ii -DBidx );
            % fprintf(1, '#%d distance=%d, shift %5.3f\n', ii, distance, unsatshifts(distance));
            chainFCH2( ii ) = chainFCH2( ii )+ unsatshifts(distance);    
        end

        % if omega < 5 the CH3 gets shifted too
        if distance <  5 
            fppmCH3 = fppmCH3 + unsatshifts(distance+1);
            % Or, from Frost & Gunstone Chem Phys Lipis 15 p53
            omegashiftsCH3 = [ 0.885, 0.890, 0.905, 0.910, 0.955];
            DBpos = sort( DBpos );

            if  (DBpos(1) <= 7) && (DBpos(1)> 2)
                fppmCH3 = omegashiftsCH3( 8-DBpos(1) );
                % fprintf( 1, 'Omega %d : %5.3f\n', DBpos(1), fppmCH3);   
            end
        end

        kdo = find( ~( DBcarbons | alphaDBcarbons) ); % excluce double bonds and the adjacent coupled CH2
        kdo = kdo(1:end-1); % exclude the last CH2. This is coupled with CH3 and included in sysFAend 

    
        %% Double bond
        %Linoleic acid
        %J = 10.7;J1=8.2;J2=-1.2;J3=6.8;J4=-2.2;

        %Oleic acid
        J = 10.6;J1=7;J2=-1.6;J3=7;J4=-1.6;

        sysFAdb.J = zeros(6);
        sysFAdb.shifts = [  5.3600
                            5.3600
                            1.9850
                            1.9850
                            1.9850
                            1.8500];

        sysFAdb.J(1, 2:6) = [J, J1, J1, J4, J4];
        sysFAdb.J(2, 3:6) = [J2, J2, J3, J3];
        sysFAdb.name = 'FAdb';
        sysFAdb.scaleFactor  = 3; % three FA per triglyceride
        % add this to the model
        sys = [sys, sysFAdb];          


        if any(alphaDBcarbons==2) 
            diallylicpos = find(alphaDBcarbons==2);
            % set the shifts for cis-cis HC=CH-CH2-HC=CH-
            shiftdiallylic = 2.665;
            % Frost and Gunstone have at 220MHz fppm = [2.720, 2.665, 2.610]
            % fhz = (fppm-4.7)*220 = -435.6000 -447.7000 -459.8000
            % Thus J = 447.7000-435.6000 or 459.8000-447.7000 = 12.1Hz!!!!!
            shiftVinyl = 5.36;
            JAX = 12.1;

            sysFAdiallylic.J = zeros(4);
            sysFAdiallylic.shifts = [   shiftdiallylic
                                        shiftdiallylic
                                        shiftVinyl
                                        shiftVinyl ];
            sysFAdiallylic.J(1, 3:4) = [JAX, JAX];
            sysFAdiallylic.J(2, 4) = JAX;

            sysFAdiallylic.name = 'FAdiallylic';
            sysFAdiallylic.scaleFactor = 3*length( diallylicpos ); % one or more diallylic, three FA per triglyceride     
            sys = [sys, sysFAdiallylic];      
        end    
    else
        % All CH2 positions switched on
        kdo = true( 1, FAchainlen-2);
        kdo = kdo(1:end-1); % exclude the last CH2. This is coupled with CH3 and included in sysFAend 
    end % length DBpos > 0
    
    % We do all the CH2, but break them up in groups of maximum 5 
    % this avoids oversized matrices in Hamiltonians
    nCH2s = length( CH2s(kdo));   
    for nn= 1:5:nCH2s
        M = nn+(0:4);
        [lia1,~] = ismember( M, CH2s);
        M = M( lia1 );
        sysFAchain.J = zeros( length(M) );
        sysFAchain.shifts = chainFCH2(M);
        sysFAchain.name = 'FAchain';
        sysFAchain.scaleFactor  = 3*2; % each two protons, three FA per triglyceride    
        sys = [sys, sysFAchain];
    end
    
    % The last CH2 is coupled with the methyl 
    sysFAend.J = zeros(5);
    sysFAend.J(1:3, 4:5) = 6.9;
    sysFAend.shifts = [   fppmCH3
                          fppmCH3
                          fppmCH3
                    chainFCH2( FAchainlen-2 )
                    chainFCH2( FAchainlen-2 )];

    sysFAend.name = 'FAend';
    sysFAend.scaleFactor  = 3; % three FA per triglyceride
    % add this to the model
    sys = [sys, sysFAend];          

    
    %% Glycerol
    J32 = 11.4; % [Hz]
    J12 = 11.9; %Hz
    J1ab = 5.88; % [Hz]
    J3ab = 5.88; % [Hz]
    s = [4.15, 4.30, 4.15, 4.30, 5.26]';
    sysGlcerol.shifts = s;

    J = zeros( 5 );
    J(1, 1:5) = [0, J1ab,   0,   0, J12];
    J(2, 1:5) = [0,   0,   0,   0, J12];
    J(3, 1:5) = [0,   0,   0, J3ab, J32];
    J(4, 1:5) = [0,   0,   0,   0, J32];

    sysGlcerol.J = J;
    sysGlcerol.name = 'Glycerol';
    sysGlcerol.scaleFactor = 1;
    % add this to the model
    sys = [sys, sysGlcerol];  
end
