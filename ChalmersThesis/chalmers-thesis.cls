% Thesis class template for Chalmers University of Technology.
% This is an inofficial class made by Mikael Öhman <mikael.ohman@chalmers.se>
% This document is released under public domain

\NeedsTeXFormat{LaTeX2e}

\ProvidesClass{chalmers-thesis}[2010/09/09 CHALMERS Thesis class]

\DeclareOption{masters}{\OptionNotUsed}
\DeclareOption{licenciate}{\OptionNotUsed}
\DeclareOption{doctorate}{\OptionNotUsed}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}

% execute default options
\ExecuteOptions{a4paper,12pt,twosided,final}

\ProcessOptions\relax

\LoadClass[a4paper,12pt,twosided]{article}

%\RequirePackage{calc}
%\RequirePackage{indentfirst}
\RequirePackage{tikz}
\RequirePackage{fancyhdr}
\RequirePackage{graphicx}

% Setting up the marginals a bit larger
\setlength{\oddsidemargin}{0mm}
\setlength{\textheight}{730pt} % default 609pt, ~710 if you dont use a4wide
\setlength{\headsep}{0pt} % default 25pt, since we have no headers in this document
\setlength{\headheight}{0pt} % 12 pt
\setlength{\voffset}{-0.4in} % default 0

\newcommand{\coverfigure}[1]{\gdef\@coverfigure{#1}}

% Introducing some new lengths for the front page
\newlength{\headerheight};
\newlength{\headeroffset};
\setlength{\headerheight}{0.215\paperheight}%~6.4cm
\setlength{\headeroffset}{0.125\headerheight}%~0.8cm

\renewcommand{\maketitle}{%
 \newpage
 \thispagestyle{empty}
 \begin{tikzpicture}[remember picture,overlay]
  \node[yshift=-6.4cm] at (current page.north west)
   {\begin{tikzpicture}[remember picture, overlay] 
     % First draw the grid and then the Logo and avancez logo.
     \clip[draw] (0,\headerheight) -- +(\paperwidth,0) -- (\paperwidth,0) -- (0.725\paperwidth,0)
                 .. controls (0.703\paperwidth,0) and (0.703\paperwidth,\headeroffset) ..
                 (0.68\paperwidth,\headeroffset) -- (0,\headeroffset) -- cycle;
     %\draw[step=5mm,black] (0,0) grid (\paperwidth,6.4cm);
     %\definecolor{logo_color}{gray}{0}
     %\node[anchor=west,xshift=0.1\paperwidth,yshift=0.511\headerheight]{\resizebox{0.63\paperwidth}{!}{\input{figures/CHALMERS.pgf}}}; % 13.2cm
     %\node[anchor=west,xshift=0.745\paperwidth,yshift=0.516\headerheight]{\includegraphics[width=0.171\paperwidth]{figures/Avancez_bg.pdf}}; % 3.6cm
   \end{tikzpicture}};
 \end{tikzpicture}
 \ \vfill
 \begin{center}
  \includegraphics[width=\textwidth,height=0.4\paperheight,keepaspectratio]{\@coverfigure}
 \end{center}
}