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
\ExecuteOptions{masters}

\ProcessOptions\relax

\LoadClass[a4paper,12pt,twosided]{article}

\RequirePackage{tikz}
\RequirePackage{fancyhdr}
\RequirePackage{hyperref}

\RequirePackage{graphicx}

% Setting up the marginals a bit larger
%\setlength{\oddsidemargin}{0mm}
%\setlength{\textheight}{730pt} % default 609pt, ~710 if you dont use a4wide
%\setlength{\headsep}{0pt} % default 25pt, since we have no headers in this document
%\setlength{\headheight}{0pt} % 12 pt
%\setlength{\voffset}{-0.4in} % default 0

%\setlength{\textwidth}{160mm}
%\setlength{\intextsep}{12pt} % This is what you meant to set to 12pt, not textfloatsep
%\setlength{\abovedisplayskip}{6pt}
%\setlength{\belowdisplayskip}{6pt}
%\setlength{\parindent}{5mm}
%\setlength{\abovecaptionskip}{6pt}
%\setlength{\belowcaptionskip}{0pt} % Actual skip = \textfloatsep + \belowcaptionskip
%\setlength{\itemsep}{0pt}

% How to nunmber different objects
%\numberwithin{equation}{section}
%\numberwithin{figure}{section}
%\numberwithin{table}{section}

\newcommand{\coverfigure}[1]{\gdef\@coverfigure{#1}}
\newcommand{\covercaption}[1]{\gdef\@covercaption{#1}}
%\newcommand{\title}[1]{\gdef\@title{#1}}
\newcommand{\subtitle}[1]{\gdef\@subtitle{#1}}
\renewcommand{\author}[1]{\gdef\@author{#1}}
\newcommand{\myabstract}[1]{\gdef\@myabstract{#1}}
\newcommand{\keywords}[1]{\gdef\@keywords{#1}}

\newcommand{\writtenyear}[1]{\gdef\@writtenyear{#1}}
\newcommand{\YYYYNN}[1]{\gdef\@YYYYNN{#1}}
\newcommand{\ISSN}[1]{\gdef\@ISSN{#1}}
\newcommand{\printers}[1]{\gdef\@printers{#1}}

\newcommand{\telephone}[1]{\gdef\@telephone{#1}}
\newcommand{\postcode}[1]{\gdef\@postcode{#1}}
\newcommand{\city}[1]{\gdef\@city{#1}}
\newcommand{\country}[1]{\gdef\@country{#1}}
\newcommand{\university}[1]{\gdef\@university{#1}}
\newcommand{\department}[1]{\gdef\@department{#1}}
\newcommand{\division}[1]{\gdef\@division{#1}}
\newcommand{\thesis}[1]{\gdef\@thesis{#1}}

% Some default settings for Chalmers
\postcode{SE-412 96}
\city{G\"oteborg}
\country{Sweden}
\university{Chalmers University of Technology}
\ISSN{1652-8557}
\printers{Chalmers Reproservice}
\telephone{+ 46 (0)31-772 1000}

% Introducing some new lengths for the front page
\newlength{\headerheight}
\newlength{\headeroffset}
\setlength{\headerheight}{0.215\paperheight}%~6.4cm
\setlength{\headeroffset}{0.125\headerheight}%~0.8cm


\AtBeginDocument{
 \hypersetup{
  breaklinks,
  pdfcreator    = \LaTeX{},
  pdfproducer   = \LaTeX{},
  bookmarksopen = true,
  pdfauthor     = \@author,
  pdftitle      = \@title,
  pdfsubject    = \@thesis,
  pdfkeywords   = \@keywords,
  %colorlinks    = true,
  urlcolor      = black,
  citecolor     = black,
  linkcolor     = black,
  filecolor     = black
 }
}

\renewcommand{\maketitle}{%
 \newpage\pagestyle{empty}
 \begin{tikzpicture}[remember picture,overlay]
  \node[yshift=-\headerheight] at (current page.north west)
   {\begin{tikzpicture}[remember picture, overlay]
     % First draw the grid and then the Logo and avancez logo.
     \clip[draw] (0,\headerheight) -- +(\paperwidth,0) -- (\paperwidth,0) -- (0.725\paperwidth,0)
                 .. controls (0.703\paperwidth,0) and (0.703\paperwidth,\headeroffset) ..
                 (0.68\paperwidth,\headeroffset) -- (0,\headeroffset) -- cycle;
     \draw[step=5mm,black] (0,0) grid (\paperwidth,6.4cm);
     \definecolor{logo_color}{gray}{0}
     \node[anchor=west,xshift=0.1\paperwidth,yshift=0.511\headerheight]{\resizebox{0.63\paperwidth}{!}{\input{figures/CHALMERS.pgf}}}; % 13.2cm
     \node[anchor=west,xshift=0.745\paperwidth,yshift=0.516\headerheight]{\includegraphics[width=0.171\paperwidth]{figures/Avancez_bg.pdf}}; % 3.6cm
   \end{tikzpicture}};
 \end{tikzpicture}
 \ \vfill
 \begin{center}
  \includegraphics[width=\textwidth,height=0.4\paperheight,keepaspectratio]{\@coverfigure}
 \end{center}
 {\huge\noindent\@title\par} % title, 21
 {\large\noindent\@subtitle\par} % subtitle, 16
 {\normalsize\noindent\textit{\@thesis}\par} % 14
 \vskip 5mm
 {\Large\noindent \uppercase\expandafter{\@author}\par}% 18
 \vskip 6mm
 {\small\noindent\@department\\
  \textit{\@division}\\
  \uppercase\expandafter{\@university}\\
  \@city, \@country\ \@writtenyear\\
  \@thesis \@YYYYNN\par} % TODO

 % Second page, blank
 \newpage
 \mbox{}

 % Third page, more info
 \newpage\thispagestyle{empty}
 \begin{center}
  {\uppercase\expandafter{\@thesis}\ \@YYYYNN\par} % WRONG
  \vspace*{\stretch{1}}%\vskip 40mm
  {\Large \@title\par}
  \vskip 5mm
  %\makeatletter\@ifundefined{subtitle}{}{{\mysubtitle\par}}\makeatother
  {\@subtitle\par}
  \vskip 5mm
  {\@thesis\par} % (if applicable) (if applicable)
  {\uppercase\expandafter{\@author}\par} % Wrong?
  \vspace*{\stretch{3}}%\vfill
  {\@department\par}
  {\textit{\@division}\par}
  {\uppercase\expandafter{\@university}\par}
  \vskip 2mm
  {\@city, \@country\ \@writtenyear\par}
 \end{center}

 % Fourth page
 \newpage
 {\noindent \@title\\
  \@subtitle
  {\uppercase\expandafter{\@author}}\par} % TODO
 \vskip 10mm
 {\noindent \copyright {\uppercase\expandafter{\@author}}, \@writtenyear\par} % TODO
 \vspace*{\stretch{1}} %\vskip 20mm
 {\noindent \@thesis\ \@YYYYNN\\ % TODO
  ISSN \@ISSN\\
  \@department\\
  \@division\\
  \@university\\
  \@postcode\ \@city\\
  \@country\\
  Telephone: \@telephone\par}
 \vspace*{\stretch{3}}
 %\makeatletter\@ifundefined{covercaption}{}{{\noindent Cover:\\\covercaption\par}\vskip 5mm}\makeatother
 {\noindent Cover:\\ \@covercaption\par} % TODO
 \vskip 5mm
 {\noindent \@printers
  %\@department
  \\ \@city, \@country\ \@writtenyear\par}

 % Fifth page
 \newpage%\pagestyle{justpage} % TODO
 \setcounter{page}{1}
 \pagenumbering{Roman}
 {\noindent \@title\\
  %\makeatletter\@ifundefined{mysubtitle}{}{\mysubtitle\\}\makeatother
  \@subtitle\\
  \@thesis\\ % (if applicable) % Thesis in programmename
  {\uppercase\expandafter{\@author}}\\
  \@department\\
  \@division\\
  \@university\par}
 \phantomsection\addcontentsline{toc}{section}{Abstract}\begin{abstract}\@myabstract\end{abstract}
 \noindent Keywords: \@keywords

}

