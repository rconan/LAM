%> @file plot_decorrelation_modes.m
%> @brief  Makes the plots for the decorrelation functions of a [modes X timeVector] covariance matrix
%> @author Mikhail Konnik
%> @date   19 May 2015
%======================================================================
%> @param correl_all     = The matrix of covariance.
%> @retval zernProj	 = The number of Zernike modes to plot
%> @retval end_point	 = how many points to plot.
% ======================================================================
function plot_decorrelation_modes(correl_all, zernProj, end_point)

number_of_modes = length(zernProj);

%%% BEGIN: Doing many plots at one
figure, plot(correl_all(1,1:end_point)); hold on;

for ii=2:number_of_modes
    plot(correl_all(ii,1:end_point));  hold on;
end

hold off
%%% END: Doing many plots at one