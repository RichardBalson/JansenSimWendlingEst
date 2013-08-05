% script created by Richard Balson 05/08/13
% This script specifies all neccesary parameters for the estimation of the
% Wendling model using simulation data from the Jansen and Rit model.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% General parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filter_simulation =0; % Specify whether or not to filter simulated data

Data_scale =1; % Scale simulation output to determine the effect of scaling on estimation out = outOrig/Data_scale

if filter_simulation
    highcutoff = 2.5; % Specify highcutoff frequency for filter
    
    lowcutoff = 40; % Specify low cutoff frequency for filter
    % Data will have frequency content between highcutoff and lowcutoff
end

fs =2048; % Specify sampling frequency

Simulation_number = 20; % Specify number of simulations to estimate

if Simulation_number>1
    Decimate = 500; % Specify the distance between corresponding samples for the output matrix when multiple simulation are performed
else Decimate =1;
end

Random_number_generator=[0 0]; % Specify whether random numbers for simulation should be different or the same for numerous simulations;

StepbyStepCheck =0; % Check how prediction and correction steps are working on a plot

LimitEst = 0; % If set to 1 estimation results are limited by physiological parameters within the estimation procedure, that is sigma points cannot be aboe the physiological parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Simulation decisions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

simulate =1; % Decide whether to simulate model output or use previous results if simulate is equal to 1 then observations for estimatio purposes are resimulated

NoiseIn = 1e-3;% Base 1e-2 Specify noise to add to simulated signal
SimulationSettings.name = 'Jansen_output'; % Specify name of file to save Wendling model output data to, or to load data from when simulation is not performed
if simulate % Specify parameters for simulation purposes
    SimulationSettings.simulation_time =100; %Time for simulation in seconds 
    SimulationSettings.slope_time =1; % Specifies the time over which the model gain should be altered
    SimulationSettings.number_of_sigma_input = 1; % Used to determine standard deviation of input if  1: 68.27% of realisations within physiolgical range, 2: 95.45, 3: 99.73 4: 99.994
    SimulationSettings.stochastic = 1; % Used to specifiy the stochastic adjustment on the input 1 is no adjustment. <1 downscalling, >1 upscaling
    SimulationSettings.Parameter_index = 4; % Choose parameters to be simulated: 1 = Seizure Parameter from Wendling 2002;
    %  2 = Seizure Parameter from Wendling 2005;...
    %  3 = Altered excitability;
    %  4 = Parameters at midpoint of their range;
    %  5 = random parameters;
    %  6= random parameters and random number of
    %  variations in simulation
    % 7 User defined For Jansen Simulation G parameter ignored
    if SimulationSettings.Parameter_index==7 % Specify synaptic gains for simulation
        SimulationSettings.AV= [5 7 5 5 5];
        SimulationSettings.BV= [20 20 20 30 20];
    end
    SimulationSettings.Input_mean_variation = 0; % If 0 mean stays constant for simulation,
    %if 1 input mean is drawn from a uniform distribution limited by the physiological limits of the input
    % if 2 input mean is drawn from a Gaussian
    % distribution with a mean as per Wendling
    % 2002 and standard deviation that satisfies
    % number_of_sigma_input
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Estimation decisions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MeanCheckTStart = 5; % For multiple simulations determine where to start checking the estimation values

EstStart = 5; % Specify the duration after simulation start when estimation should start. This allows removal of all transients.

% Estimation states
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ds = 8; % Number of differential equations describing model, also the number of fast states to be estiamted

Dp = 3; % Number of parameters to be estimated, also refered to as slow states

Dk =1; %If set to 1 the mean of the stochastic input will be estimated % Note that if Input_mean_variation is not zero than Dk should be set to one to allow tracking of the input mean

Dy =1; % Number of observable outputs from the simulation

% Intialisation of parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Parameter_initialisation = 0; % 0 for parameters to be initialised by Gauss distribution
% 1 for parameters to be initialsed as a
% percetage error from actual values

if Parameter_initialisation
    PercError = 10; % Specify percentage error for parameter intialisation
end

number_of_sigma = 1; % Number of standard deviations from mean. 4 accounts for 99.73 percent of points.

Reinitialise_parameters_attempts = 1; % Specify number of attempts for parameter reinitialisation if results are not physiologucally possible

% Uncertainty parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

kappa =0; % Varibale used to define the relative contribution of the mean on the propogation of states, and adjustment of the variance of the sigma points drawn from the Gaussian distribution

Variable_state_uncertainty = 0;%1e-3; % 1e-3 Uncertianty due to stochastic input

State_uncertainty_adjustment = [1 2 3 40 1 2 3 40];% Exponential decrease in uncertainty % All ones good for slow but steady convergence

Exc_parameter_uncertainty = 1e-2;
SInh_parameter_uncertainty =8e-2;
FInh_parameter_uncertainty =8e-2;
Base_parameter_uncertainty = 8e-2;%1e-2;%1e-6;%1e-3; % Inherent parameter uncertainty due to model error

Variable_parameter_uncertainty = 0;%1e-3;  % Uncertianty due parameters varying in time

Base_input_uncertainty = 1e-12;%1e-3; % Inherent parameter uncertainty due to model error

Variable_input_uncertainty =0;%1e-3; % Uncertianty due varying input mean, Set to zero if the input mean is not varying

Observation_uncertainty = 1e-3; % Specify the uncertainty in observations

uncertainty_adjustment = 1; % Adjuster for model uncertainty

std_adjustment_parameters =1; %2% Variance adjuster for parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Image handling parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Estimation_Type = ['GF',int2str(filter_simulation),'_PS',int2str(SimulationSettings.slope_time)]; % Estimation Type is an indication of what estimation is being performed for, here Gauss indicates that all staes are initialised as realisations from a Gaussian distribution
                           % Also used for the saving of figures % GF Gauss
                           % filter, PS parameter slope
                           
fig_save =1; % Save figures as .fig for future use

printpdf =1; % Specify whether results should be printed to pdf

Image_handling_model_output=[1;0];

plot_uncertainty =1; % Plot covariance of all states

Image_handling_states = [0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0]; % Here a decision is made whether to plot specific states, 
                                                                            % if the value is one the relevant figure is plotted, otherwise it is not.
                                                                            % The columns indicate the state to plot and the rows indicate whether the whole simulation or a zoomed in ploted should be plotted.
                                                                            % Column one corresponds with state 1 and so forth.


Image_handling_multi = [1 1;0 0];%  % Here a decision is made whether to plot specific states, 
                                                                            % if the value is one the relevant figure is plotted, otherwise it is not.
                                                                            % The columns indicate the figures to be plotted.
                                                                            % Here column 1-2 are for all the model states and all the model parameters including the input mean.

plot_uncertaintyMulti =0;
% Zoom parameters (seconds)
% ~~~~~~~~~~~~~~~~~

tstart =0; % Starting time for zoom

zoomtime = 10; % Duration of zoom