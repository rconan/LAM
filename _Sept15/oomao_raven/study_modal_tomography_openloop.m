reconstruction_method = 'AmmseLag';


%% Begin Open Loop
tel = tel+atm;
coefsVec = zeros(dm.nValidActuator,number_of_methods2copmare,nIteration);
nSlope = wfs(1).nSlope;

for kIteration = 1:nIteration
    fprintf('Iteration %d out of %d \n',kIteration,nIteration)
    
    for iFrame = 1:frameTime/samplingTime

        +tel;  % Always advance atm by 1ms

        for iGs = 1:nGs;
            ast(iGs) = ast(iGs).*tel*wfs(iGs);
        end

            sciObjectVec{1} = sciObjectVec{1}.*tel*dmObjectVec{1}; %%% ???
            sciObjectVec{1} = sciObjectVec{1}*imgrObjectVec{1};


        if iFrame == fixedLagTimeInMs
            if kIteration == 1
                    dmObjectVec{1}.coefs = -coefsVec(:,1,kIteration);
            else
                    dmObjectVec{1}.coefs = -coefsVec(:,1,kIteration-1);
            end
        end

    end %%% END generating the frames


    for iGs = 1:nGs;
        slopesStack(:,iGs,kIteration) = wfs(iGs).slopes;
    end

         NGS(1:nSlope,kIteration) = slopesStack(:,1,kIteration);

%%% This is for 3 NGSes
%         NGS(1:nSlope,kIteration) = slopesStack(:,1,kIteration);
%         NGS(nSlope+1:nSlope*2,kIteration) = slopesStack(:,3,kIteration);
%         NGS(2*nSlope+1:nSlope*3,kIteration) = slopesStack(:,2,kIteration);


    %% Reconstruction methods
    switch reconstruction_method

         case 'Static MMSE'
%%%% Static MMSE reconstruction
             coefsVec(:,1,kIteration) = Z2U * Rmv * NGS(:,kIteration);     %%% The problem appears to be here: the slopes are a lot larger than the coefficientsVectors

         case 'AR1'
%%%% AR1 prediction
            coefsVec(:,1,kIteration) = Z2U*RmvAR1*NGS(:,kIteration);

         case 'AmmseLag'
            coefsVec(:,1,kIteration) = Z2U*RmvAmmse*NGS(:,kIteration);

    end

end %%% for kIteration

%% Why there is no instantaneous Strehl ratio available?
imgrObjectVec{1}.ee     %% encircled energy
imgrObjectVec{1}.strehl %% strehl ration
% ./OOMAO-Raven/OOMAOlibUpdated/imager.m