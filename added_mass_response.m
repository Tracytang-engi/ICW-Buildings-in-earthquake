%% Added mass location in a three-storey structure (NO absorber)
% Adapted from the Cambridge Moodle three-storey MATLAB model.
% Added mass is rigidly attached to a floor: the system always has 3 DOFs.
% Input: 1 N harmonic horizontal force on floor 1; base fixed.
% c = 3.96 Ns/m references the A1 Python code's primary damping DEFAULT.
% Applying that value to every interstorey connection is our assumption;
% it is not a measured damping calibration of this three-storey model.
clear; close all; clc;

%% 1. Parameters and output folder
m = 1.83;                 % Original mass of EACH floor, kg
dm = 0.5*m;               % Representative added mass: 0.915 kg
c = 3.96;                 % EACH connection's viscous damping, Ns/m
L = 0.2; b = 0.08; E = 210e9; d = 0.001;
I = b*d^3/12;
k = 24*E*I/L^3;           % 4200 N/m
fmin = 0.1; fmax = 20;    % Hz, not rad/s
df = 0.0025;              % Display grid; maxima are refined below
f = linspace(fmin,fmax,round((fmax-fmin)/df)+1);
F = [1;0;0];              % Force amplitude in N
S = [2 -1 0; -1 2 -1; 0 -1 1];
K = k*S;
C = c*S;                 % Includes base-to-floor-1 damping
caseNames = {'Baseline';'Added floor 1';'Added floor 2';'Added floor 3'};
caseShort = {'Baseline','Floor 1','Floor 2','Floor 3'};
col = lines(4);
styles = {'-','--','-.',':'};
outdir = fullfile(pwd,['added_mass_results_' datestr(now,'yyyymmdd_HHMMSS')]);
if ~exist(outdir,'dir'), mkdir(outdir); end

%% 2. Mass matrices, undamped frequencies and mass-normalised modes
Mcase = zeros(3,3,4);
Phi = zeros(3,3,4);
fn = zeros(3,4);
zeta = zeros(3,4);
eigResidual = zeros(4,1);
for s = 1:4
    masses = m*ones(3,1);
    if s > 1, masses(s-1) = masses(s-1) + dm; end
    Mcase(:,:,s) = diag(masses);
    [fn(:,s),Phi(:,:,s)] = getModes(K,Mcase(:,:,s));
    wn = 2*pi*fn(:,s);
    zeta(:,s) = c*wn/(2*k);
    eigResidual(s) = norm(K*Phi(:,:,s) - ...
        Mcase(:,:,s)*Phi(:,:,s)*diag(wn.^2),'fro') / ...
        norm(K*Phi(:,:,s),'fro');
end
freqChange = 100*(fn./fn(:,1)-1);
frequencyTable = table(caseNames,fn(1,:)',fn(2,:)',fn(3,:)', ...
    freqChange(1,:)',freqChange(2,:)',freqChange(3,:)', ...
    'VariableNames',{'Case','Mode1_Hz','Mode2_Hz','Mode3_Hz', ...
    'Mode1_change_pct','Mode2_change_pct','Mode3_change_pct'});
disp('Undamped natural frequencies and changes:'); disp(frequencyTable);
writetable(frequencyTable,fullfile(outdir,'natural_frequencies.csv'));
dampingTable = table(caseNames,100*zeta(1,:)',100*zeta(2,:)', ...
    100*zeta(3,:)','VariableNames', ...
    {'Case','Mode1_zeta_pct','Mode2_zeta_pct','Mode3_zeta_pct'});
writetable(dampingTable,fullfile(outdir,'modal_damping_ratios.csv'));

%% 3. LOCAL mass sensitivity at the equal-mass baseline
% Modes are mass normalised: phi_j' * M * phi_j = 1.
% df_j/dm_i = -f_j*phi_ij^2/2, in Hz/kg.
% This is a derivative, not the exact shift for dm = 0.5*m.
dfdm = zeros(3,3);         % Rows: modes; columns: added-mass floors
for j = 1:3
    dfdm(j,:) = -0.5*fn(j,1)*(Phi(:,j,1).^2)';
end
relativeSensitivity = (m./fn(:,1)).*dfdm;
modeNames = {'Mode1';'Mode2';'Mode3'};
Ts = array2table(dfdm,'VariableNames',{'Floor1','Floor2','Floor3'}, ...
    'RowNames',modeNames);
writetable(Ts,fullfile(outdir,'sensitivity_Hz_per_kg.csv'),'WriteRowNames',true);
Tr = array2table(relativeSensitivity, ...
    'VariableNames',{'Floor1','Floor2','Floor3'},'RowNames',modeNames);
writetable(Tr,fullfile(outdir,'sensitivity_dimensionless.csv'),'WriteRowNames',true);
% Independent central-difference check of the sensitivity formula.
h = 1e-4*m;
fdSensitivity = zeros(3,3);
for i = 1:3
    Mp = Mcase(:,:,1); Mm = Mp;
    Mp(i,i) = Mp(i,i)+h; Mm(i,i) = Mm(i,i)-h;
    fp = getModes(K,Mp); fm = getModes(K,Mm);
    fdSensitivity(:,i) = (fp-fm)/(2*h);
end
sensitivityError = max(abs(fdSensitivity-dfdm),[],'all') / ...
    max(abs(dfdm),[],'all');

%% 4. Damped harmonic responses: direct dynamic-stiffness solution
X = complex(zeros(3,numel(f),4));
peakByFloor = zeros(3,4);
peakHzByFloor = zeros(3,4);
modalAgreement = zeros(4,1);
gridGap = zeros(4,1);
for s = 1:4
    M = Mcase(:,:,s);
    for q = 1:numel(f)
        w = 2*pi*f(q);
        X(:,q,s) = (K-w^2*M+1i*w*C)\F;
    end
    % Cross-check direct solution against modal superposition.
    Xm = modalResponse(f,fn(:,s),Phi(:,:,s),c/k,F);
    modalAgreement(s) = max(abs(X(:,:,s)-Xm),[],'all') / ...
        max(abs(X(:,:,s)),[],'all');
    % Refine each observed local maximum, so a peak between grid points
    % is not mistaken for a reduction in response.
    [peakByFloor(:,s),peakHzByFloor(:,s)] = refinePeaks( ...
        f,abs(X(:,:,s)),M,C,K,F);
    gridGap(s) = max((peakByFloor(:,s)-max(abs(X(:,:,s)),[],2)) ...
        ./peakByFloor(:,s));
end
A = abs(X);              % Absolute displacement amplitudes, metres
[worstPeak,criticalFloor] = max(peakByFloor,[],1);
criticalHz = zeros(1,4);
for s = 1:4
    criticalHz(s) = peakHzByFloor(criticalFloor(s),s);
end
peakChange = 100*(worstPeak/worstPeak(1)-1); % Positive means WORSE
peakTable = table(caseNames,1e3*peakByFloor(1,:)', ...
    1e3*peakByFloor(2,:)',1e3*peakByFloor(3,:)', ...
    1e3*worstPeak',criticalFloor',criticalHz',peakChange', ...
    'VariableNames',{'Case','Floor1_peak_mm','Floor2_peak_mm', ...
    'Floor3_peak_mm','Worst_peak_mm','Critical_floor', ...
    'Critical_frequency_Hz','Worst_peak_change_pct'});
disp('Positive peak change means a larger maximum response:'); disp(peakTable);
writetable(peakTable,fullfile(outdir,'peak_responses.csv'));
% Export the full sampled response curves with phase preserved in the MAT file.
frfTable = table(f','VariableNames',{'Frequency_Hz'});
for s = 1:4
    for i = 1:3
        name = sprintf('Case%d_Floor%d_amplitude_m',s-1,i);
        frfTable.(name) = reshape(A(i,:,s),[],1);
    end
end
writetable(frfTable,fullfile(outdir,'frequency_response_amplitudes.csv'));

%% 5. Sweep added mass; retain the restriction that no mass is removed
% 0-100%: main range; 110-500%: broader numerical sensitivity check.
% 500% is an exploratory model calculation, not a physical load rating.
massRatios = unique([0:0.01:1, 1.1:0.1:5])';
sweepPeak = zeros(numel(massRatios),3);
sweepCriticalHz = sweepPeak;
fCoarse = linspace(fmin,fmax,1991);
for i = 1:3
    for r = 1:numel(massRatios)
        masses = m*ones(3,1);
        masses(i) = masses(i)+massRatios(r)*m;
        M = diag(masses);
        [fLocal,P] = getModes(K,M);
        % Extra samples around every natural frequency resolve narrow peaks.
        nearModes = fLocal.*linspace(0.9,1.1,81);
        fs = unique([fCoarse,nearModes(:)']);
        fs = fs(fs>=fmin & fs<=fmax);
        Xs = modalResponse(fs,fLocal,P,c/k,F);
        [pk,pf] = refinePeaks(fs,abs(Xs),M,C,K,F);
        [sweepPeak(r,i),imax] = max(pk);
        sweepCriticalHz(r,i) = pf(imax);
    end
    fprintf('Added-mass sweep completed for floor %d.\n',i);
end
sweepChange = 100*(sweepPeak/worstPeak(1)-1);
sweepTable = table(massRatios,m*massRatios, ...
    1e3*sweepPeak(:,1),1e3*sweepPeak(:,2),1e3*sweepPeak(:,3), ...
    sweepChange(:,1),sweepChange(:,2),sweepChange(:,3), ...
    'VariableNames',{'Added_mass_ratio','Added_mass_kg', ...
    'Added_floor1_peak_mm','Added_floor2_peak_mm','Added_floor3_peak_mm', ...
    'Added_floor1_change_pct','Added_floor2_change_pct','Added_floor3_change_pct'});
writetable(sweepTable,fullfile(outdir,'added_mass_sweep.csv'));
[bestPeak,indexBest] = min(sweepPeak,[],1);
bestTable = table((1:3)',m*massRatios(indexBest(:)),1e3*bestPeak', ...
    100*(bestPeak'/worstPeak(1)-1),'VariableNames', ...
    {'Added_floor','Best_sampled_added_mass_kg','Worst_peak_mm','Change_pct'});
disp('Best sampled masses, including the zero-added-mass option:'); disp(bestTable);
writetable(bestTable,fullfile(outdir,'best_sampled_mass.csv'));

%% 6. Figures: export both PNG and editable MATLAB FIG
% Four separate undamped frequency-response plots
floorColours = lines(3);

for s = 1:4
    M = Mcase(:,:,s);
    Xu = nan(3,numel(f));

    for q = 1:numel(f)
        % Skip points extremely close to an undamped natural frequency
        if min(abs(f(q)-fn(:,s))) < 1e-6
            continue
        end

        w = 2*pi*f(q);
        Xu(:,q) = (K-w^2*M)\F;
    end

    % Break the plotted curves at each singularity
    for j = 1:3
        [~,idx] = min(abs(f-fn(j,s)));
        Xu(:,idx) = NaN;
    end

    g1 = figure('Color','w','Position',[80 80 1000 550]);
    hold on;

    for i = 1:3
        plot(f,1e3*abs(Xu(i,:)), ...
            'Color',floorColours(i,:), ...
            'LineWidth',1.4, ...
            'DisplayName',sprintf('Response of floor %d',i));
    end

    set(gca,'YScale','log');
    grid on;
    xlim([fmin fmax]);

    xlabel('Excitation frequency (Hz)');
    ylabel('Displacement amplitude (mm)');
    title(['Undamped frequency response: ' caseNames{s}]);
    subtitle('Read resonance frequencies; amplitudes are unbounded at resonance');
    legend('Location','northeast');

    saveFigure(g1,outdir,sprintf('01_undamped_case%d',s-1));
end

g2 = figure('Color','w','Position',[80 80 1100 420]);
tiledlayout(1,3,'TileSpacing','compact');
for j = 1:3
    nexttile; hold on;
    for s = 1:4
        v = Phi(:,j,s)/max(abs(Phi(:,j,s)));
        % Align signs with the baseline; eigenvector sign is arbitrary.
        if Phi(:,j,1)'*v < 0, v = -v; end
        plot([0;v],(0:3)',styles{s},'Color',col(s,:), ...
            'LineWidth',1.5,'Marker','o');
    end
    xline(0,':'); grid on; yticks(0:3); ylim([0 3]); xlim([-1.15 1.15]);
    xlabel('Normalised modal displacement'); ylabel('Floor (0 = fixed base)');
    title(sprintf('Mode %d',j));
    if j==1, legend(caseNames,'Location','best'); end
end
saveFigure(g2,outdir,'02_mode_shapes');

g3 = figure('Color','w','Position',[80 40 1050 850]);
tiledlayout(3,1,'TileSpacing','compact');
for i = 1:3
    nexttile; hold on;
    for s = 1:4
        plot(f,1e3*reshape(A(i,:,s),1,[]),styles{s}, ...
            'Color',col(s,:),'LineWidth',1.4);
    end
    set(gca,'YScale','log');
    grid on; xlim([fmin fmax]); ylabel('Amplitude (mm, log scale)');
    title(sprintf('Response of floor %d to a 1 N force on floor 1',i));
    if i==1, legend(caseNames,'Location','northeast'); end
    if i==3, xlabel('Excitation frequency (Hz)'); end
end
saveFigure(g3,outdir,'03_damped_frequency_responses');

g4 = figure('Color','w','Position',[80 80 1050 430]);
tiledlayout(1,2,'TileSpacing','compact');
for panel = 1:2
    nexttile; hold on;
    for i = 1:3
        plot(100*massRatios,sweepChange(:,i),styles{i+1}, ...
            'Color',col(i+1,:),'LineWidth',1.6);
    end
    yline(0,'k:'); grid on;
    xlabel('Added mass / original single-floor mass (%)');
    ylabel('Change in band-wide maximum displacement (%)');
    if panel==1
        xlim([0 100]); title('Main added-mass range');
        legend(caseNames(2:4),'Location','northwest');
    else
        xlim([0 500]); title('Extended numerical check');
    end
end
saveFigure(g4,outdir,'04_added_mass_sweep');

%% 7. Verification and reproducibility record
staticX = K\F;
fprintf('\nStatic displacements (mm): %.6f %.6f %.6f\n',1e3*staticX);
fprintf('Maximum relative eigen-equation residual: %.3g\n',max(eigResidual));
fprintf('Relative sensitivity check error: %.3g\n',sensitivityError);
fprintf('Maximum direct/modal response disagreement: %.3g\n',max(modalAgreement));
fprintf('Largest sampled/refined peak difference: %.4f %%\n',100*max(gridGap));
if max(eigResidual)>1e-8 || sensitivityError>1e-5 || max(modalAgreement)>1e-8
    warning('A numerical cross-check exceeded its tolerance. Inspect results.');
end
save(fullfile(outdir,'all_results.mat'), ...
    'm','dm','c','k','L','b','E','d','I','F','K','C','Mcase', ...
    'f','fmin','fmax','df','fn','Phi','zeta','dfdm','relativeSensitivity', ...
    'X','A','peakByFloor','peakHzByFloor','worstPeak','peakChange', ...
    'massRatios','sweepPeak','sweepChange','sweepCriticalHz', ...
    'eigResidual','sensitivityError','modalAgreement','gridGap', ...
    'frequencyTable','dampingTable','peakTable','sweepTable','bestTable');
thisFile = [mfilename('fullpath') '.m'];
if exist(thisFile,'file'), copyfile(thisFile,fullfile(outdir,'script_used.m')); end
fid = fopen(fullfile(outdir,'run_notes.txt'),'w');
fprintf(fid,'MATLAB: %s\nDate: %s\n',version,datestr(now));
fprintf(fid,'m = %.8g kg; dm = %.8g kg; c = %.8g Ns/m; k = %.8g N/m\n',m,dm,c,k);
fprintf(fid,'Frequency band %.8g to %.8g Hz; force amplitude 1 N on floor 1.\n',fmin,fmax);
fprintf(fid,'No absorber. Fixed base. Added mass is rigidly attached.\n');
fprintf(fid,'c is a modelling assumption referencing the A1 default, not measured here.\n');
fprintf(fid,'Positive peak change means worse response. Sweep minima are sampled minima.\n');
fprintf(fid,'Eigen residual %.4g; sensitivity error %.4g; modal agreement %.4g.\n', ...
    max(eigResidual),sensitivityError,max(modalAgreement));
fclose(fid);
fprintf('\nAll figures, tables, code and raw results saved in:\n%s\n',outdir);

%% Local functions (keep these at the end of the file)
function [fn,P] = getModes(K,M)
    [P,D] = eig(K,M);
    [lambda,order] = sort(real(diag(D)),'ascend');
    P = real(P(:,order));
    fn = sqrt(lambda)/(2*pi);
    for j = 1:size(P,2)
        P(:,j) = P(:,j)/sqrt(P(:,j)'*M*P(:,j));
        [~,imax] = max(abs(P(:,j)));
        if P(imax,j)<0, P(:,j) = -P(:,j); end
    end
end

function X = modalResponse(f,fn,P,beta,F)
    wn = 2*pi*fn(:); w = 2*pi*f(:)';
    denominator = wn.^2-w.^2+1i*beta*(wn.^2)*w;
    X = P*((P'*F)./denominator);
end

function [pk,pf] = refinePeaks(f,A,M,C,K,F)
    [pk,idx] = max(A,[],2);
    pf = reshape(f(idx),[],1);
    opts = optimset('TolX',1e-9,'Display','off');
    for i = 1:3
        candidates = find(A(i,2:end-1)>=A(i,1:end-2) & ...
            A(i,2:end-1)>=A(i,3:end))+1;
        for q = candidates
            [fq,negativeAmplitude] = fminbnd( ...
                @(ff)-floorAmplitude(ff,i,M,C,K,F),f(q-1),f(q+1),opts);
            if -negativeAmplitude>pk(i)
                pk(i) = -negativeAmplitude; pf(i) = fq;
            end
        end
    end
end

function a = floorAmplitude(f,i,M,C,K,F)
    w = 2*pi*f;
    x = (K-w^2*M+1i*w*C)\F;
    a = abs(x(i));
end

function saveFigure(g,outdir,name)
    exportgraphics(g,fullfile(outdir,[name '.png']),'Resolution',300);
    savefig(g,fullfile(outdir,[name '.fig']));
end
