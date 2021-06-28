Fig.\ \ref{curves} shows the viral titer curves for different MOI of SARS-CoV-2, where the darker line for each color shows the curve of median values and the lighter colored lines are the 100 individual simulations. Note that for most MOI, there is very little variation between simulations once the viral titer is large. The exception is the lowest MOI of 10$^{-5}$ where there is more variation in the exact trajectory of the viral load. We see some obvious shifts in the viral titer curve as the MOI increases. For high MOI, the viral titer curve reaches its peak very quickly, with lower MOIs moving the peak farther out in time. The peak also becomes broader and lower as the MOI becomes lower, suggesting longer infection durations, but with lower viral loads.
\begin{figure}[!h]
\begin{center}
\resizebox{0.6\textwidth}{!}{\includegraphics{Figures/Covid_AllOnOne_Median_VirusVsTime}}
\caption{Viral loads for infections initiated with different MOI. Dark lines of each color indicate the viral load curve using the median of 100 simulations, while the lighter colored lines show the viral load kinetics for each individual simulation. The dashed line indicates the threshold of detection used to calculate infection duration. \label{curves}}
\end{center}
\end{figure}

For a more quantitative assessment, we measure the characteristics described in Methods. The results are shown in Fig.\ \ref{results}, which shows peak viral load (top left), time of viral peak (top right), viral upslope (center left), AUC (center right), and infection duration (bottom) as functions of the MOI. The peak viral load increases with increasing initial inoculum, but it appears to reach a plateau as we near an MOI of 1. The time of peak, on the other hand, decreases with increasing initial inoculum, reaching a fixed small value at high MOI. There are real plateaus here since each cell will produce an average of $p\tau_I$ viral particles. At an MOI of 1, all cells are initially infected and will start producing virus at about the same time, meaning that all that virus is released almost simultaneously and there is no second cycle of infection. At slightly lower MOIs, most cells are initially infected, but some cells will be infected in a second or third cycle of infection, reducing the large burst of virus at one time, which causes a delay, reduction, and broadening in the peak. The upslope, or growth rate, of the viral titer curve increases as the MOI increases. This is also driven by the larger amount of virus being produced in the first cycle of infection as the MOI increases. Finally, the AUC and infection duration both decrease as the initial inoculum increases.  
\begin{figure}[!h]
\begin{center}
%\resizebox{0.48\textwidth}{!}{\includegraphics{Figures/CovidApectGraphs/PeakViralTitter}}
%\resizebox{0.48\textwidth}{!}{\includegraphics{Figures/CovidApectGraphs/TimeofPeakViralTitter}}
%\resizebox{0.48\textwidth}{!}{\includegraphics{Figures/CovidApectGraphs/New_UpSlope}}
%\resizebox{0.48\textwidth}{!}{\includegraphics{Figures/CovidApectGraphs/AUC}}
%\resizebox{0.48\textwidth}{!}{\includegraphics{Figures/CovidApectGraphs/InfectionDuration}}
\caption{Effect of initial inoculum on viral titer characteristics. The graphs show peak viral load (top left), time of viral peak (top right), viral upslope (center left), AUC (center right), and infection duration (bottom) as functions of MOI. \label{results}}
\end{center}
\end{figure}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Discussion}

Our study finds that initial viral inoculum does alter the viral time course by increasing the peak viral load, moving the peak earlier, increasing the viral upslope, and decreasing both AUC and infection duration, as the initial inoculum increases. It is not immediately clear what these changes in viral kinetics mean for the severity of the infection. Is it better to have a shorter infection, albeit with a higher viral peak, or a longer-lasting infection with a lower viral burden? One study compared viral loads in patients with mild and severe illness finding that the viral load time course in mild cases peaked earlier and at a lower peak viral load than in severe cases, although both time courses still had rather high viral loads at 25 days post symptom onset \cite{zheng20}. Since viral load in these patients was measured after they presented at a hospital, there is also no way to link particular features of the viral time course to the initial inoculum. Other observational studies that have attempted to investigate links between viral load and disease severity have taken a limited number of viral load measurements, often well after the peak of the infection \cite{liu20, liu20imm,to20}, making it impossible to assess the full time course of the viral load, and any correlations to initial inoculum. An alternative to observational studies in patients is to investigate inoculum dose-response of SARS-CoV-2 in animals, as suggested in \cite{little20}. Such animal studies in conjunction with mathematical modeling studies will help provide a clearer picture of the role of initial inoculum in determining viral time course and disease severity.
