%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Viruses are microscopic parasites, generally much smaller than bacteria. They lack the capacity to thrive and reproduce outside of a host body \cite{website2}. A virus is composed of a nucleic acid genome and a protein capsid that covers the genome \cite{website3}. As seen in figure \ref{fig:Virus_Replication}, the life cycle of a virus begins with the virus attaching to or being absorbed by the host cell. Once the virus genome enters into a cell, the genome moves to the ribosomes, where the genome is replicated. With the replicated genome, new virus can be assembled and released from the host cell \cite{Kaiser}, allowing for them to continue spreading through out the host cells.

\begin{figure}[h]
    \centering
    \includegraphics[width=0.6\linewidth]{Figures/Virus_Replication_large.pdf}
    \caption{The life cycle of a virus begins with a virion (virus particle) being absorbed by the cell. Once the virion enters into a cell the virus genome is released. The genome moves to the ribosomes, where the genome is replicated. With the replicated genome, new virus can be assembled by the golgi apparatus and then released from the cell.}
    \label{fig:Virus_Replication}
\end{figure}

Some viruses cause illnesses and a few of them are severe enough to receive global recognition: the 2019-2020 coronavirus pandemic (a widespread global outbreak), the 2014 outbreak of Ebola in West Africa, and the 2009 H1N1/swine flu pandemic, for example. Another virus that has become well known is influenza (or the flu). In total, the Centers for Disease Control and Prevention (CDC) estimates that in the United States up to 42.9 million people were sick during the 2018-2019 flu season, 647,000 people were hospitalized, and 61,200 died \cite{website4}.

In order to understand viruses, assays are performed. An assay is an experiment for assessing or measuring characteristics of a substance. The two typical forms of assays are quantitative and quantal. Quantitative assays are assays that give an accurate and exact measure of the amount of a substance in a sample. Of these types of assays, plaque assays are the most widely used for determining viral titer \cite{Kumar}. They can be used with any virus that causes damage to the cells where they have been grown. This damage is called a plaque and is circular in shape. These assays are often performed in petri dishes, where virus is placed in a dish of healthy cells and the formation of plaques and the concentration of virus are monitored. It is assumed that each plaque formed is caused by one virus particle. Because of this assumption, the viral concentration is often recorded as plaque forming units per milliliter (pfu/mL). Quantal assays are assays which generally give only a pass or fail. These assays are performed in many mediums such as in animals or in multi-well plates. Any virus that will infect the medium can be used. When multi-well plates are used, different groups of wells on the plate are filled with a dilution of virus that often varies in a ten fold difference from group to group. After a set amount of time the wells are observed for a reaction, and whether or not a well has a response is counted. When animals are used, different animal subjects are injected with a dilution of virus that often varies in a ten fold difference from subject to subject. After a set amount of time the animals are observed, and whether or not an animal dies is counted.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Agent-based (individual-based or micro-simulation) models have been around since 1970 with the introduction of ``Conway's Game of Life'' \cite{gardner70}. These models have been utilized in many different fields from physics to the study of fish (ichthyology) \cite{owusu20} and continue to be popularized for different applications by software like Netlogo \cite{nogare20,chiacchio14}. The models consist of a collection of agents whose behavior is determined by mathematical or computational rules. The agents of the system can move freely \cite{beauchemin07} or be fixed in a grid or lattice \cite{beauchemin_simple_2005} for varying applications, but either configuration allows for tracking of spatially emergent patterns. In recent years, the field of virology has started using agent-based models to study the spread of viruses in a monolayer of cells \cite{beauchemin_simple_2005,alvarado_cellular-level_2018,wodarz_laws_2014,tong_development_2015,whitman20,goyal16,itakura10,wasik14} in an effort to study the spatiality of viral spread.

%In a lab, the most common monolayer of cells is a petri dish. Petri dishes are a type of adherent culture where the cells are grown on a nutritious substrate. A petri dish is circular in shape and often has a diameter between 35-150 mm. When the cells are grown to confluence the cells tend to push on each other and distort the shape of each cell membrane \cite{bruckner_importance_2018}.
In a lab, \emph{in vitro} viral infections are performed on layers of cells grown to the point of confluence, where there is on the order of \numrange[range-phrase = --]{e5}{e6} cells \cite{Number_of_cells_in_a_dish_noauthor_useful_nodate}. Virus modelers are using ABMs to simulate the two dimensional layer of agents to replicate experiments that are done \emph{in vitro} in order to better understand factors that affect viral spread. The agent-based model framework is appealing to virus modelers because it allows for the individual tracking of how cells, as agents, interact with the virus, and has the potential to replicate \emph{in vitro} and eventually \emph{in vivo} viral infections. Currently, however, the implementation of agent-based models in the field of virology has two issues: speed and size. 

Agent-based models are notorious for being computationally intensive and taking long amounts of time to run simulations. This point has been commented on in a previous article \cite{gallagher_causes_2018} and the feasibility of ABMs for research has been talked about as a goal that is to come with increasing computational advancements \cite{bauer_agent-based_2009}. Previous research has addressed this lack of computing power issue by reducing the number of agents modeled and therefore reducing the number of computations required for a simulation. The number of agents published is at minimum an order of magnitude lower then the number of target cells used in the corresponding experimental data. Beauchemin et al.\ \cite{beauchemin_simple_2005} simulated \num{1.232e5} agents, while the experiment they were attempting to replicate was performed in 6 well-plates and had $\sim$\num{1.2e6} cells per well. Alvarado et al.\ \cite{alvarado_cellular-level_2018} simulated \num{4.0e4} agents when trying to replicate experiments also performed in 6 well-plates. Wodarz et al.\ \cite{wodarz_laws_2014} simulated \num{2.0e4} agents, while the experiment they were replicating was performed in 24 well-plates and had $\sim$\num{2.4e5} cells per well. Tong et al.\ \cite{tong_development_2015} simulated \num{6.0e5} agents in an effort to simulate mice lungs, which have $\sim$\num{e9} cells. These smaller simulations are more affected by boundary interactions, which can result in model dynamics that don't faithfully reproduce the infection. This can hinder the research process and can make it difficult to use models for treatment optimization and personalized medicine.

While it might be feasible to wait long periods of time to run accurate simulations for endemic or recurrent seasonal viruses, recent events of the COVID-19 pandemic indicate how great a need there is for accurate and fast modeling methods. Epidemiological population-level modeling tools that include both ordinary differential equation (ODE) models \cite{li20,ngonghala20} and ABMs \cite{ying21,sneppen21,kano21} were immediately deployed to help predict how the new virus would spread around the world and how different interventions could help stem the spread. At the within-host level, the primary modeling tool was limited to simple ODE models \cite{goncalves20,wang20model,hernandez20,dogra20} that lack the ability to reproduce the spatial heterogeneity of real viral infections. Fast and accurate in-host models could be helpful in assessing the potential of re-purposed drugs \cite{czuppon21,goncalves20,dodds20}, finding indicators of disease severity or mortality \cite{neant21}, and assessing the effectiveness of testing \cite{ejima21}. A community-driven ABM incorporating many realistic biological responses was quickly developed for SARS-CoV-2 \cite{getz21}, but is only currently simulating a few thousand agents and is expected to need high-performance computing or cloud resources to parameterize the model. Thus, there is a need to develop modeling and simulation tools for accurately predicting in-host viral dynamics that can be quickly deployed to help combat the next pandemic.

In this paper, we present the testing and validation of a hybrid ABM and partial differential equation model (PDEM) implemented on GPUs, which allows for rapid simulation of large-scale simulations on desktop computers. The work here begins with the methods where the four attributes of the model, the agent-based model of the cells, the partial differential equation of the virus, the cell-free transmission mode of viruses, and fitting of the model to data, are explained. We then present the results of model implementation with parallel programming, convergence testing, and simulation speed improvement. Finally, we show that the model can reproduce experiments by fitting the model to an example data set from an \emph{in vitro} influenza experiment.

\subsection{Model details} 

In this work, a two dimensional biological system is simulated with a mathematical model. The system is a culture dish of a monolayer of cells with virus diffusing over the cells. The model is a hybrid of an agent based model (ABM) and a partial differential equation model (PDEM) where the cells are represented with an ABM and the virus diffusion is represented by a PDEM.

\subsubsection{Viral transmission} \label{Viral_transmission}

When a virus is spreading among the cells in a culture dish, there is a probability that a healthy cell becomes infected by virus that is not within a cell, but flowing around and above the cell. When this viral transmission occurs it is called cell-free transmission. For cell-free transmission, the probability per unit time ($\mathrm{P_{cf}}$) that a cell becomes infected is determined by the amount of virus that is covering the cell ($V$) times the infection rate ($\beta$) \cite{holder11autoimm}, 
$$\mathrm{P_{cf}} = V \beta.$$ 
As a healthy cell becomes surrounded by more virus, the probability of cell-free infection increases. If the probability ($V \beta \Delta t$) is even greater than one due to the build up of virus, an adaptive time step is used. The time step ($\Delta t$) is divided in half repeatedly until the probability of cell-free infection is below one. Once the probability is finalized, a number from the uniform distribution is compared with the probability of cell-free infection. If that number is less than $\mathrm{P_{cf}}$, then the cell becomes infected.

%As the viral infection progresses, the total amount of virus in the culture dish changes. Plotting the amount of virus vs.\ time produces a curve that has a distinct shape and has characteristics that can be measured. The measurements, shown in Figure \ref{measurements}, can be used to compare multiple viruses or to compare multiple simulations of the same virus with different input parameters.

%\begin{itemize}
%\item \textbf{peak viral load:} The maximum amount of virus is commonly used as an indicator of the transmissibility of an infection \citep{handel09}. 
%\item \textbf{time of viral peak:} This is the time between the start of the infection and the peak of the virus and can give an indication of how quickly the virus is replicating.
%\item \textbf{viral upslope:} Viral upslope is the exponential growth rate of the viral titer before the peak is reached and is another indication of how quickly the virus is spreading from cell to cell. 
%\item \textbf{viral downslope:} \color{red}Text missing here.\color{black}
%\item \textbf{area under the curve (AUC):} AUC is often used to assess the severity of an infection \citep{hayden00, barroso05}.
%\item \textbf{infection duration:} The infection duration is indicative of how long an infected patient might test positive for presence of the virus. In this work $10^1 \mathrm{PFU/ml}$ is the threshold.
%\end{itemize}
%\begin{figure}
%\begin{center}
%    \resizebox{0.6\textwidth}{!}{\includegraphics{Figures/measurements.pdf}}
%    \caption{Measurable characteristics of the viral titer curve.\color{red} Figure is missing downslope label.\color{black}}
%    \label{measurements}
%\end{center}
%\end{figure}

%1 - y = (1-m)^n

\subsubsection{Spatial accounting} \label{Spatial_accounting}

To allow for the two dimensional aspect of the culture dish to be represented in the model, the cells are approximated as hexagons. Using hexagons enables for an elegant managing of the cells' shapes in the dish and the viral transmission. Since the culture dishes are grown to confluence, the cells are close enough that they push on each other and the cell walls deform. This causes the cells to no longer be in the shape of a circle, but become irregular polygons with multiple sides \cite{bruckner_importance_2018}. Modeling the cells as hexagons gives the cells definite sides and the cells are able to span any two dimensional region forming a hexagonal grid. Furthermore, by using a hexagonal grid, when virus particles spread among this population of cells the indexing of the grid can be used to find the neighbors of any cell. This will be used for cell-free transmission to know where virus will flow away from (high concentrations areas) and to (low concentration areas) during diffusion. 

In addition to helping with the physical representation of the model, hexagonal coordinates have some other attributes that can be utilized to optimize the code for quicker compute times. The three attributes that this code utilizes are:
\begin{enumerate} 
    \item The coordinates can be split in to three sectors where the coordinates $X_{hex}$, $Y_{hex}$, and $Z_{hex}$ are simply cyclic permutations.
    \item The $X_{hex}$ and $Z_{hex}$ directions can be used as indices of a matrix.
    \item The coordinates of the neighboring hexagons are found by adding a cyclic permutation of 
        $\left [
            \begin{array}{c}
                1 \\
                0 \\
                -1\\
            \end{array}
        \right ]$
        for three of the neighbors and
        $\left [ 
            \begin{array}{c}
                1 \\
                -1 \\
                0\\
            \end{array}
        \right ]$
        for the other three neighbors.
\end{enumerate}
These attributes save time by either reducing the number of calculations needed or the amount of searching through data arrays. Attribute 1 allows for only a third of the cell locations to be calculated and Attributes 2 and 3 give the data a reference so that adjacent data in memory can be found quicker.

\subsubsection{Agent-based model} \label{ABM}

In an ABM, a system is broken down into smaller units called ``agents''. Each of the agents are governed by a set of rules on a local scale with large scale phenomena resulting from interaction of the agents, so the two scales are studied to find the connections. As a simulation of the model is stepped through time, the agents act and interact. These actions cause bulk properties, that may appear disconnected from the individual agents, to manifest. Properties are observed and measured to find the connection between the small interactions and large scale properties.

In this work, an ABM governs the transitions a cell makes through the stages of infection: healthy, eclipse, infected, and dead. A cell in the healthy state is an uninfected cell that remains healthy until infected. A cell in the eclipse state is an infected cell that is not yet producing virus. The cell remains in the eclipse state for an average amount of time, $\tau_E$. The specific time value for each cell is determined by a gamma distribution with shape value $\eta_E$ and scale value $\tau_E/\eta_E$. A cell in the infected state is an infected cell that is producing virus. The cell remains in the infected state for an average amount of time, $\tau_I$. The specific time value for each cell is determined by a gamma distribution with shape value $\eta_I$ and scale value $\tau_I/\eta_I$. A gamma (Erlang) distribution is used for the amount of time in the eclipse and infected phase, as suggested by the work of Beauchemin et al.\ \cite{beauchemin17} and Kakizoe et al.\ \cite{kakizoe15}. A cell in the dead state is a cell that can no longer change state, so once a cell is in the dead state the cell remains in that state until the end of the simulation. The flow of this is illustrated in figure \ref{fig:transitioning_through_the_stages_of_infection}.

The ABM uses four time arrays to track and transition the cells to different states after infection. The four arrays universal time (UT), healthy time (HT), eclipse time (ET), and infected time (IT) have an element for each cell. The universal time array holds the amount of time that each cell has been in the simulation; each element starts at zero and increases each iteration by the simulation's time step. The healthy time array holds the amount of time that a cell is healthy; each element starts at zero and while the cell is healthy increases each iteration by the simulation's time step. The eclipse time array holds the amount of time each cell is in the eclipse state and the infected time array holds the amount of time each cell is in the infected state. For the eclipse and infected arrays the amount of time is fixed and the value is determined by a gamma (Erlang) distribution, as described above. The flow of this is illustrated in figure \ref{fig:transitioning_through_the_stages_of_infection}.

%\begin{figure}
%    \centering
%    \begin{tikzpicture}[state/.style={regular polygon,regular polygon sides=6, draw, minimum size=2.5cm, inner sep=0pt, outer sep=0pt}]
%    \node[state, draw=AGreen, fill=AGreen!10] (H) {$Healthy$};
%    \node[state, draw=cyan, fill=cyan!10, right of=H] (E) {$Eclipse$};
%    \node[state, draw=red, fill=red!10, right of=E] (I) {$Infected$};
%    \node[state, draw=black, fill=black!10, accepting, right of=I] (D) {$Dead$};
%    
%    \draw   (H) edge[bend left=45] node [above] {} (E);
%    \draw   (E) edge[bend left=45] node [above] {$\tau_E$} (I);
%    \draw   (I) edge[bend left=45] node [above] {$\tau_I$} (D);
%    \end{tikzpicture}
%    \caption{The stages of infection healthy, eclipse, infected, and dead are shown. Cells stay in  the eclipse stage for an average time $\tau_E$ and in the infected stage for an average time $\tau_I$.}
%    \label{fig:stages_of_infection}
%\end{figure}

%\begin{figure}
%    \centering
%    \begin{tikzpicture}[state/.style={regular polygon,regular polygon sides=6, draw, minimum size=2.5cm, inner sep=0pt, outer sep=0pt}]
%        \node[state, draw=AGreen, fill=AGreen!10] (H) {$Healthy$};
%        \node[state, draw=cyan, fill=cyan!10, below right of=H] (E) {$Eclipse$};
%        \node[state, draw=red, fill=red!10, below right of=E] (I) {$Infected$};
%        \node[state, draw=black, fill=black!10, accepting, below right of=I] (D) {$Dead$};
%            
%        \draw   (H) edge[bend left=45] node [above right] {Infection event} (E);
%        \draw   (E) edge[bend left=45] node [above right] {UT $>$ HT + ET} (I);
%        \draw   (I) edge[bend left=45] node [above right] {UT $>$ HT + ET + IT} (D);
%    \end{tikzpicture}
%    \caption{Transitioning through the stages of infection: UT is the universal time for a cell, HT is the healthy time for a cell, ET is the eclipse time for a cell, and IT is the infected time for a cell.}
%    \label{fig:transitioning_through_the_stages_of_infection}
%\end{figure}

\begin{figure}
    \centering
    \begin{tikzpicture}[
        Healthy/.style={regular polygon,regular polygon sides=6, draw, minimum size=2.5cm, inner sep=0pt, outer sep=0pt, draw=AGreen, fill=AGreen!10},
        Eclipse/.style={regular polygon,regular polygon sides=6, draw, minimum size=2.5cm, inner sep=0pt, outer sep=0pt, draw=cyan, fill=cyan!10},
        Infected/.style={regular polygon,regular polygon sides=6, draw, minimum size=2.5cm, inner sep=0pt, outer sep=0pt, draw=red, fill=red!10},
        Dead/.style={regular polygon,regular polygon sides=6, draw, minimum size=2.5cm, inner sep=0pt, outer sep=0pt, draw=black, fill=black!10}]

        \node[Healthy] (H) {$Healthy$};
        \node[Eclipse, below right of=H] (E) {$Eclipse$};
        \node[Infected, below right of=E] (I) {$Infected$};
        \node[Dead, accepting, below right of=I] (D) {$Dead$};
            
        \draw   (H) edge[bend left=45] node [above right] {Infection event} (E);
        \draw   (E) edge[bend left=45] node [above right] {UT $>$ HT + ET} (I);
        \draw   (I) edge[bend left=45] node [above right] {UT $>$ HT + ET + IT} (D);

        \draw   (H) edge[bend right=45] node [below left] {} (E);
        \draw   (E) edge[bend right=45] node [below left] {$\tau_E$} (I);
        \draw   (I) edge[bend right=45] node [below left] {$\tau_I$} (D);

        %Legand with text
        \matrix [draw,
                below right,
                column 1/.style={anchor=base west},
                column 2/.style={anchor=base west}
                ] at (current bounding box.north east) 
        {
          \node [] {UT}; & \node [] {Universal time}; \\ 
          \node [] {HT}; & \node [] {Healthy time}; \\ 
          \node [] {ET}; & \node [] {Eclipse time}; \\ 
          \node [] {IT}; & \node [] {Infected time}; \\ 
        };

%        %Legand with nodes
%        \matrix [draw, below left] at (current bounding box.north east) {
%          \node [Healthy, scale=0.3, label=right:Healthy] {}; \\
%          \node [Eclipse, scale=0.3, label=right:Eclipse] {}; \\
%          \node [Infected, scale=0.3, label=right:Infected] {}; \\
%          \node [Dead, scale=0.3, label=right:Dead] {}; \\
%        };
    \end{tikzpicture}
\caption{The stages of infection: healthy, eclipse, infected, and dead are shown. The cells transition through the stages at different time points. Above: The time point when a state transition occurs is shown in terms of UT, the universal time, for a cell. %UT is the time a cell has existed, HT is the time a cell has been healthy, ET is the time a cell is in the eclipse phase, and IT is the time a cell is in the  infected.
Below: The time point when a state transition occurs is shown in terms of average time. $\tau_E$ is the average time a cell stays in  the eclipse stage and $\tau_I$ is the average time a cell stays in the infected stage. \label{fig:transitioning_through_the_stages_of_infection}}
\end{figure}

\subsubsection{Partial differential equation model} \label{PDEM}

PDEMs are used to model multiple dimensions; in this work a PDE in hexagonal coordinates is used to model the two-dimensional spatial spread of virus over cells in a culture dish. In a PDEM, the dynamics of a system can be represented by a partial differential equation, or more specifically, an equation that contains multi-variable functions that represent important system aspects and one or more partial derivatives of those functions. In the culture dish, as an infected cell releases virus into the extracellular fluid, the virus diffuses across a density gradient. The PDEM represents this diffusion with the diffusion equation, 
\begin{equation}
\frac{\partial V}{\partial t}=D \nabla^{2}V + p - cV, \label{diff_eq}
\end{equation}
where $V$ is the density of the virus, $D$ the diffusion coefficient, $p$ is the production rate per cell, $c$ is the viral clearance rate. In the code, along with the assumption of hexagonal cells, the cells are assumed to be flat, so the virus is diffusing over a smooth two dimensional plane. This assumption allows for the use of the two dimensional diffusion equation in hexagonal coordinates, so Eq.\ \eqref{diff_eq} becomes 
$$\frac{\partial V}{\partial t} = D\frac{2}{3} \left (\frac{\partial^2}{\partial x^2_1}+\frac{\partial^2}{\partial x^2_2}+\frac{\partial^2}{\partial x^2_3}\right )V + p -cV$$ 
where  
$\textbf{x}_1=
\left [
    \begin{array}{c}
        1 \\
        0 \\
    \end{array}
\right ]$, 
$\textbf{x}_2=
\left [
    \begin{array}{c}
        -1/2 \\
        \sqrt{3}/2 \\
    \end{array}
\right ]$, and 
$\textbf{x}_3=
\left [
    \begin{array}{c}
        -1/2 \\
        -\sqrt{3}/2 \\
    \end{array}
\right ]$ 
are the unit vectors for a hexagonal grid. For computation, we use a forward Euler implementation of the PDEM with Neumann boundary conditions.

\subsubsection{Parameters of viral spread}

The eight parameters $\beta$, $\tau_E$, $\eta_E$, $\tau_I$, $\eta_I$, $p$, $c$, and $D$ affect the dynamics of virus spread in the model. Four of the parameters, $\tau_E$, $\eta_E$, $\tau_I$, and $\eta_I$, are used in the ABM to choose the time duration that a cell is in the eclipse and infected phase as mentioned in section \ref{ABM}. Three of the other parameters, $p$, $c$, and $D$, are used in the PDEM and characterize the differential equation, as mentioned in section \ref{PDEM}. The final parameter, $\beta$, governs the interaction between the virus and cells, setting the probability that the cell is infected. In order to model a particular virus, values for these parameters need to be chosen. The initial values of the parameters are chosen from ordinary differential equation models of influenza and listed in Table \ref{tab_params} (viral titer units have been converted to virions, as described in \cite{dobrovolny17}). %For $\eta_E$, $\eta_I$, $c$, and $D$ the values are fixed, but for $\beta$, $\tau_E$, $\tau_I$, and $p$ the values serve as a starting point for the model to be fit to real experimental data.

\begin{table}
%\centering
\caption{Parameter values to simulate an influenza infection with the ABM/PDEM model.\label{tab_params}}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llcr}
\hline
Parameter & Meaning & Value & Reference\\
\hline
$\beta$ & Infection rate & 2.0 $/\mathrm{h}$ & Scaled from Beauchemin et al.\ \cite{beauchemin08}\\
$p$ & Viral production rate & 562800 $/\mathrm{h}$ & Scaled from Beauchemin et al.\ \cite{beauchemin08}\\
$c$ & Viral clearance rate & 0.105 $/\mathrm{h}$ & Beauchemin et al.\ \cite{beauchemin08}\\
$D$ & Diffusion coefficient & 2.16$\times 10^{-8}$ $\mathrm{m}^2/\mathrm{h}$ & Stokes-Einstein equation\\
$\tau_E$ & Mean eclipse duration & 6.0 $\mathrm{h}$ & Beauchemin et al.\ \cite{beauchemin08}\\
$\eta_E$ & Eclipse shape parameter & 30 & Pinilla et al.\ \cite{pinilla12}\\
$\tau_I$ & Mean infectious lifespan & 12.0 $\mathrm{h}$ & Beauchemin et al.\ \cite{beauchemin08}\\
$\eta_I$ & Infectious shape parameter & 100 & Pinilla et al.\ \cite{pinilla12}\\
\end{tabular}}
\end{table}

%In order to determine which parameter values minimize the SSR, a brute force walk of a parameter space is performed. The axes of the parameter space are four arrays that are constructed around the initial values of $\beta$, $p$, $\tau_I$, and $\tau_E$, shown in Table \ref{tab_params}. Two values, increasing or decreasing by 25\%, are added on both sides of the initial points; resulting in a total of 5 points per array and 625 points total in the parameter space.

\subsection{Computational details}

\subsubsection{Implementation on GPUs}

As the model becomes more complex, GPU acceleration via parallel programming is used to decrease the simulation run times and therefore increase the number of studies that can be conducted in a given time. In the simulations, the cells change state based on the amount of virus above them. The number of cells in a culture dish is on the order of $10^6$ cells \cite{Number_of_cells_in_a_dish_noauthor_useful_nodate}, so the ABM will simulate a grid of $1001365$ agents of hexagonal cells in a circle to best replicate what is happening in the experiment. Each agent will follow the rules of checking the amount of virus above the cell every time step. Utilizing attribute 2 of hexagonal coordinates, the number of calculations is reduced from the order of ($\mathcal{O}(n^2)$) per time step to the order of the number of agents ($\mathcal{O}(n)$). The calculations from the agents' rules are split over the processing units of a GPU to be calculated in parallel or simultaneously. To utilize this processing, Nvidia's CUDA (Compute Unified Device Architecture) is used to implement the ABM and PDEM. CUDA is an Application Programming Interface (API) that allows the many processing units (cores) on a Nvidia brand GPU to be used for computing. %The ability to manipulate what is happening at the cellular level of the simulations allows for generation of isolated studies of the viral transmission. The isolated studies of how cells are infected are analyzed and the characteristics from section \ref{Viral_transmission} are compare to find any trends in the data.

\subsubsection{Convergence Testing}

%Partial differential equations (PDEs) are a popular way to model systems that evolve over both space and time, an example is discussed in section \ref{PDEM}. With PDEs, even systems that have an exact solution need to be calculate on a computer, because of the infinite series that are required in those solutions. So solutions to PDEs are often found through numerical integration. In the numerical integration, space and time are assumed to be made up of small units (discretization). Time is a one dimensional line of points separated by a chunk of time called $\Delta t$ and space is a grid with a line of points for each dimension where there is a chunk of space for each dimension $\Delta x$, $\Delta y$ (two dimensional Cartesian space). At these points in time and space, a numerical integration scheme (numerical scheme) for approximating the solution of the PDE is chosen. Different numerical schemes have different benefits. Depending on the phenomena that needs to be studied with the PDE the size of $\Delta t$, $\Delta x$, and $\Delta y$ and the choice of numerical scheme are important. If the chunks of space or time are too large then the simulation does not have the resolution to resolve phenomena that occur at smaller increments in the model and if the numerical scheme requires to much computing power then the solutions can not be found in a timely manner.

%From the choice of numerical scheme a conditional relationship between $\Delta t$, $\Delta x$, and $\Delta y$ is formed and must be met. For the Euler's method $$ \frac{\Delta t}{(\Delta x)^{2}} \leq \frac{1}{2},$$ if we make $\Delta x = \Delta y$.\cite{cite} This relationship is necessary to ensure that the sequence of approximations, that the numerical scheme uses to approximate a solution, converges and not allow the error to grow exponentially to a point that the solutions are unreliable. 

%Using the relationship above, a value for $\Delta t$, $\Delta x$, and $\Delta y$ can be chosen to ensure stability of the error in the numerical scheme. As long as that relationship is met the solution is reliable within a certain error, but the relationship does not give the $\Delta t$, $\Delta x$, and $\Delta y$ that are best for producing accurate simulations with the least amount of computing cost.

To ensure the simulation of the PDEM converges, we need to optimize the space and time discretizations: $\Delta t$, $\Delta x$, and $\Delta y$. Convergence testing is a simple brute force method where the input parameters are increased or decreased by a particular amount and the accuracy or trends of the simulation are measured for each of the the new increments. Schemes for convergence testing are implemented and studied in fields like computational fluid dynamics \cite{bermejo16,kim20fluid} and astrophysics \cite{xu21,banei21}. The model proposed in this work has fixed $\Delta x$ and $\Delta y$, because the simulations are of real cells, whose average diameter can be measured between \numrange[range-phrase = --]{50}{100}\si{\micro\meter}. Thus the convergence testing only has to be conducted for $\Delta t$. To conduct the study a starting point of 0.005 hr was chosen and a range of seven values was created by multiplying or dividing the initial $\Delta t$ by 2 repeatedly. For each of these $\Delta t$s, the median viral titer curve of ten simulations were compared.

\subsection{Data Fitting} \label{Data_Fitting}

As part of our model validation, we verified that the model could reproduce viral titer curves observed experimentally. The experimental data set used here is from an \emph{in vitro} experiment performed by Pinilla et al.\ \cite{pinilla12}. During the study, a well of a 24-well plate, containing Madin-Darby canine kidney (MDCK$\alpha$2,6) cells was inoculated with the A/Qu\'{e}bec/144147/09 (H1N1) pandemic strain of influenza virus and the supernatant fluid was collected every 6 hours until 36 hours and then every 12 hours until 72 hours post infection. The supernatant was then used for RNA isolation and/or viral titration by standard plaque assay on MDCK$\alpha$2,6 cells. The specific data referenced for this work is the multiple-cycle viral yield experiment shown in figure 2A of the Pinilla manuscript.

To determine the best fit of the model to the experimental data, the sum of square residuals (SSR) is minimized, $$\mathrm{SSR} = \sum_{i=1}^{n} (y_i - \hat y_i)^{2},$$ where $y_i$ is from the experimental data set and $\hat y_i$ is from the simulated data set. In our case, the simulated data set is the average of ten cell-free transmission simulations. The initial conditions for the simulations are: Total cells -- $1001365$, Total virus -- $0.0$, and MOI -- $5\times 10^{-5}$. To perform the minimization, a separate code that utilizes the function \texttt{minimize} from the python package \texttt{scipy}, was written. In the code, five parameters ($\beta$, $p$, $\tau_I$, $\tau_E$, and $c$) are allowed to vary and the remaining parameters are held fixed to the values given in Table \ref{tab_params}. The minimization code is given an initial guess for the five parameters, then by the Nelder-Mead method the next set of parameters is produced, until the minimum SSR is found.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The novel coronavirus, SARS-CoV-2, originated in Wuhan, China in late 2019 and rapidly spread around the world \cite{chen20,wu20}. While the disease can lead to severe illness needing long hospitalization \cite{sun20,goyal20,jiang20}, a significant fraction of those who contract the virus are asymptomatic \cite{he20}. It is still not entirely clear who is at risk for developing severe disease, although correlations of disease severity with levels of vitamin D \cite{ilie20}, levels of various immune components \cite{liu20imm,liu20imm2,zhang20imm,yang20imm}, and age \cite{borghesi20,zhang20imm} have been noted. There has also been investigation of the possibility of disease severity being linked to initial viral inoculum \cite{little20, guallar20, ghandi20}.

There is some evidence from other respiratory viruses that the size of the initial inoculum could play a role in the severity of the illness. An influenza epidemiological modeling study suggests that a higher initial dose can lead to a higher mortality rate \cite{paulo10}. This is corroborated by an influenza in-host modeling study that also finds a correlation between the initial viral dose and survival rate \cite{price15}. Other modeling studies have found dependence of other measures of infection severity on initial dose for influenza \cite{moore20}, respiratory syncytial virus \cite{wethington19}, adenovirus \cite{li14}, and porcine reproductive and respiratory virus \cite{go19}. There are also experimental studies that find a link between dose and infection severity. Experiments using influenza have found inoculum dose dependence of total number of infected cells and area under the curve \cite{manicassamy10}, peak viral titer \cite{ginsberg52,iida63,ottolini05}, viral growth rate \cite{ginsberg52}, and time of viral peak \cite{iida63,ginsberg52}. Experiments with other viruses, such as adenovirus \cite{prince93}, and parainfluenza \cite{ottolini96}, have also shown correlations between initial inoculum and various measures of disease severity. If SARS-CoV-2 shows a similar pattern, initial inoculum should be considered as a possible contributor to infection severity and adverse outcomes.

The major route of transmission for SARS-CoV-2 is airborne droplets \cite{morawska20}. One study indicates that sneezing and coughing creates a turbulent gas cloud that can cause viral-laden droplets to spread up to 27 feet (\numrange[range-phrase = --]{7}{8}\si{\meter}) \cite{bourouiba20}, and allows the virus to get into the ventilation system of a building. A review of literature on droplet and airborne viral spread concludes that 8 of 10 studies showed that droplets spread further than the 6 foot \cite{bahl20} social distancing recommendation. While personal protective equipment is helpful in reducing the ability of virus to enter the respiratory tract, it is not perfect \cite{mittal20}. All of these factors lead to exposures to vastly different quantities of virus when people are going about their daily activities. Thus it is important to understand whether different initial inocula lead to different viral dynamics in patients. 

Given the difficulty of examining SARS-CoV-2 inoculum dependence in patients, our study aims to address the question of inoculum-dependence of SARS-CoV-2 infection severity using mathematical modeling. We use a combination agent-based model (ABM) and partial differential equation model (PDM) to simulate SARS-CoV-2 infections initiated with different initial inocula. We measure several features of the viral titer curve and find that increasing the initial inoculum leads to an early, high, and narrow peak in the viral titer curve, while decreasing both the infection duration and area under the curve.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\subsection{Mathematical model}

We use an ABM to model transitions of cells as they go through the infection cycle. We use a hexagonal grid and simulate 10$^6$ cells in a circular dish to mimic an in vitro system. Cells begin as healthy target cells that can be infected by viruses that are sitting above them. Once infected, the cells move into an eclipse phase where they are not yet actively producing virus. The cells remain in the eclipse phase for a time chosen from an Erlang distribution with mean time $\tau_E$ and shape parameter $n_E$. The cells then pass into the infectious phase, where they are actively producing virus, for a time chosen from an Erlang distribution with mean time $\tau_I$ and shape $n_I$, after which time the cells die and no longer participate in the infection. Erlang distributions are used for both transitions based on experiments that show the time spent in the eclipse phase and the time spent in the infectious phase are best described by Erlang distributions \cite{kakizoe15, beauchemin17}, at least for SHIV. While SHIV is a different virus, it is the only virus for which these distributions have been measured directly. Influenza, another respiratory virus, has also been shown to need non-exponential transition distributions \cite{holder11autoimm, holder11}. 

Viral dynamics are described by the PDM as virus diffuses over the layer of cells,
\begin{equation}
\frac{\partial V}{\partial t} = D\nabla^2V+p-cV,
\end{equation}
where $D$ is the diffusion coefficient, and $c$ is the viral decay rate. Virus is produced by infectious cells at rate $p$ and is assumed to be released directly above each infected cell. The amount of virus above any cell determines the probability that the cell will be infected, $P_\mathrm{inf}=\beta V$, where $P_\mathrm{inf}$ is the probability per unit time, and $\beta$ is the infection rate. A more detailed description of the model is given in the supplementary material and the simulation code is available on https://github.com/BaylorFain/Covid19-Code.

Parameter values that describe SARS-CoV-2 are taken from a variety of sources and are given in Table \ref{params}. The majority of the parameters are taken from \cite{pinky20}, where an ordinary differential equation model of coronavirus infection was fit to viral titer data from a single patient. Note that the parameters $\beta$ and $p$ are scaled to account for the different numbers of cells (10$^6$ here and 1 in \cite{pinky20}) in the two systems as well as converting viral concentration to individual virions (see \cite{handel07,perelson12,dobrovolny17} for detailed discussions on converting from concentration to virions). The shape parameters are based on values derived from influenza infections \cite{pinilla12}, since the Erlang distribution has not yet been used for SARS-CoV-2. The diffusion coefficient was calculated using the Stokes-Einstein equation \cite{cush97}. 
\begin{table}
\centering
\caption{Parameter values to simulate SARS-CoV-2 infection with the ABM/PDM model.\label{params}}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llc}
\hline
Parameter & Meaning & Value \\
\hline
$\beta^a$ & Infection rate & $\SI{84.0}{\per\hour}$ \\
$\tau_E^b$ & Mean eclipse duration & \SI{5.88}{\hour} \\
$n_E^c$ & Eclipse shape parameter & 30 \\
$\tau_I^b$ & Mean infectious lifespan & \SI{0.624}{\hour} \\
$n_I^c$ & Infectious shape parameter & 100 \\
$p^a$ & Viral production rate & $\SI{19900}{\per\hour}$ \\
$c^b$ & Viral clearance rate & \SI{0.00490}{\per\hour} \\
$D^d$ & Diffusion coefficient & $\SI{4.80e-12}{\meter^2\per\second}$ \\
\hline
\multicolumn{2}{l}{$^a$Parameters taken from \cite{pinky20}, but scaled.}\\
\multicolumn{2}{l}{$^b$Parameters taken from \cite{pinky20}.}\\
\multicolumn{2}{l}{$^c$Parameters taken from \cite{pinilla12}.}\\
\multicolumn{2}{l}{$^d$Parameter calculated from Stokes-Einstein equation.}
\end{tabular}}
\end{table}

\subsection{Measurements}

We simulate SARS-CoV-2 infections starting with different multiplicity of infection (MOI) where the MOI value defines the initial number of infected cells. The ABM/PDM model is implemented in Compute Unified Device Architecture (CUDA) and run on NVIDIA graphics processing units. We perform 100 simulated infections for each MOI and measure the following features of the viral titer curve (Fig.\ \ref{measurements}): 
\begin{itemize}
\item \textbf{peak viral load:} The maximum amount of virus is commonly used as an indicator of the transmissibility of an infection \citep{handel09}. 
\item \textbf{time of viral peak:} This is the time between the start of the infection and the peak of the virus and can give an indication of how quickly the virus is replicating.
\item \textbf{viral upslope:} Viral upslope is the exponential growth rate of the viral titer before the peak is reached and is another indication of how quickly the virus is spreading from cell to cell. 
\item \textbf{area under the curve (AUC):} AUC is often used to assess the severity of an infection \citep{hayden00, barroso05}.
\item \textbf{infection duration:} The infection duration is indicative of how long an infected patient might test positive for presence of the virus. Note that the threshold used here is 10$^7$ virions based on a 10$^2$ RNA copies/ml detection threshold for the experimental data \cite{goncalves20} that is converted to individual virions.
\end{itemize}
\begin{figure}[!h]

\begin{center}
\resizebox{0.6\textwidth}{!}{\includegraphics{Figures/measurements}}
\caption{Characteristics of the viral titer curve that are used to assess severity of the infection.\label{measurements}}
\end{center}
\end{figure}
