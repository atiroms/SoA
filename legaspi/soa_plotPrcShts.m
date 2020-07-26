% Published: August 14, 2019
% Copyright
% Lab for Neural Computation and Adaptation
% RIKEN Center for Brain Science
% Objective: Graph the action and outcome perceptual shifts as functions of
% the strenth of the causal prior

function F = soa_plotPrcShts(experiment, arrPrcShftA, arrPrcShft0, arrPXil, fontsize)

F = figure;
linewidth = 2;

if experiment == 1
    % Haggard et al., 2002 (Nat Neurosci): Seminal intentional binding experiment
    % Different keypress (i.e., the action) conditions

    plot(arrPXil,arrPrcShftA(1,:),'b--', arrPXi1,arrPrcShft0(1,:),'b--',arrPXil,arrPrcShftA(2,:),'r--', arrPXi1,arrPrcShft0(2,:),'r--',arrPXi1,arrPrcShftA(3,:),'k--', arrPXi1,arrPrcShft0(3,:),'k--', 'Linewidth', Linewidth);

    %legend('Voluntary ActionÅf, ' Tone', ÅeInvoluntary MEPÅf, ' Tone',' Sham TMS',' Tone', ÅeLocation', 'northwest');
    lgnd = legend('Voluntary action',' and tone', 'Involuntary action',' and tone','Sham',' and tone', 'Location', 'northwest', 'Orientation','vertical');
    lgnd.FontSize = 18;
    set(lgnd.BoxFace, 'ColorType', 'truecoloralpha', 'ColorData', uint8(255*[1; 1; 1; 0.8]));

elseif experiment == 2
    % Wolpe et al. 2013 (Exp Brain Res): Uncertainty is with the outcome
    % Different tone (i.e., the outcome) conditions
    hold on;
    plot(arrPXil,arrPrcShftA(1,:),'Color', [0 0 250/255], 'LineStyle','-', 'Linewidth', linewidth);
    plot(arrPXil,arrPrcShft0(1,:),'Color', [0 0 250/255], 'LineStyle', '--', 'Linewidth', linewidth) ;
    plot(arrPXil,arrPrcShftA(2,:),'Color', [0 140/255 255/255], 'LineStyle','-', 'Linewidth', Linewidth);
    plot(arrPXil,arrPrcShft0(2,:),'Color', [0 140/255 255/255], 'LineStyle', '--', 'Linewidth', Linewidth);
    plot(arrPXil,arrPrcShftA(3,:),'Color', [0 240/255 255/255], 'LineStyle','-', 'Linewidth', linewidth);
    plot(arrPXil,arrPrcShft0(3,:),'Color', [0 240/255 255/255], 'LineStyle', '--', 'Linewidth', linewidth) ;
    lgnd = legend('Action', ' and low uncertainty tone', 'Action', ' and intermediate uncertainty tone', 'Action', ' and high uncertainty tone','Location', 'northwest');
    lgnd.FontSize = 18;
    set(lgnd.BoxFace, 'ColorType', 'truecoloralpha', 'ColorData', uint8(255*[1; 1; 1; 0.8]));
end

set(gca,'FontSize', fontsize);
set(gca, 'Box', 'Åeon');
set(lgnd, 'Color', 'none');
set(gca,'color', 'white');
set(gca,'FontSize', fontsize);
