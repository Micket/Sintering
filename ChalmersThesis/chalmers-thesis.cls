% Thesis class template for Chalmers University of Technology.
% This is an inofficial class made by Mikael Öhman <mikael.ohman@chalmers.se>
% This document is released under public domain
% Some words are translated to swedish if babel langauge is set, otherwise default to English

\NeedsTeXFormat{LaTeX2e}

\ProvidesClass{chalmers-thesis}[2010/09/09 CHALMERS Thesis class]

%\DeclareOption{bachelor}{\@bachelortrue}
%\DeclareOption{masters}{\@masterstrue}
%\DeclareOption{licenciate}{\@licenciatetrue}
%\DeclareOption{doctorate}{\@doctoratetrue}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}

% execute default options
\ExecuteOptions{masters}

\ProcessOptions\relax

\LoadClass[a4paper,12pt]{article}

\RequirePackage{tikz}
\RequirePackage{fancyhdr}
\RequirePackage{hyperref}
\RequirePackage{graphicx}
%\RequirePackage{xcolor}

\definecolor{phd_thesis_blue}{rgb}{0.4, 0.4, 0.8}
\definecolor{chalmers_blue}{rgb}{0, 0, 0.4}
\definecolor{chalmers_grey}{rgb}{0.8, 0.8, 0.8}


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

%\renewcommand{\author}[1]{\gdef\@author{#1}}
%\newcommand{\title}[1]{\gdef\@title{#1}}
\newcommand{\subtitle}[1]{\gdef\@subtitle{#1}}
\newcommand{\myabstract}[1]{\gdef\@myabstract{#1}}
\newcommand{\keywords}[1]{\gdef\@keywords{#1}}
\newcommand{\coverfigure}[1]{\gdef\@coverfigure{#1}}
\newcommand{\covercaption}[1]{\gdef\@covercaption{#1}}

\newcommand{\writtenyear}[1]{\gdef\@writtenyear{#1}}
\newcommand{\YYYYNN}[1]{\gdef\@YYYYNN{#1}}
\newcommand{\ISSN}[1]{\gdef\@ISSN{#1}}
\newcommand{\ISBN}[1]{\gdef\@ISBN{#1}}
\newcommand{\printers}[1]{\gdef\@printers{#1}}

\newcommand{\telephone}[1]{\gdef\@telephone{#1}}
\newcommand{\postcode}[1]{\gdef\@postcode{#1}}
\newcommand{\city}[1]{\gdef\@city{#1}}
\newcommand{\country}[1]{\gdef\@country{#1}}
\newcommand{\university}[1]{\gdef\@university{#1}}
\newcommand{\department}[1]{\gdef\@department{#1}}
\newcommand{\division}[1]{\gdef\@division{#1}}
\newcommand{\thesis}[1]{\gdef\@thesis{#1}}

% Some default settings for Chalmers or at least the department of solid mechanics.
\postcode{SE-412 96}
\city{G\"oteborg}
%\iflanguage{swedish}{\country{Sverige}}{\country{Sweden}}
\country{\iflanguage{swedish}{Sverige}{Sweden}}
%\university{\iflanguage{swedish}{Chalmers Tekniska Högskola}{Chalmers University of Technology}}
\university{Chalmers University of Technology} % expandafter just won't work, MakeUppercase is annoying.
\printers{Chalmers Reproservice}
\telephone{+ 46 (0)31-772 1000}
%\if@bachelor
  \ISSN{1654-4676} % Bachelor
%\fi
\ISSN{1652-8557} % Masters
\ISSN{1652-8565} % Lic.
\ISSN{0346-718X} % PhD

\newif\if@restonecol
\newif\if@titlepage
\@titlepagefalse


% Introducing some new lengths for the front page
\newlength{\headerheight}
\newlength{\headeroffset}
\setlength{\headerheight}{0.215\paperheight}%~6.4cm
\setlength{\headeroffset}{0.125\headerheight}%~0.8cm


\AtBeginDocument{
{
 \renewcommand{\and}{, }
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
}

\renewcommand{\maketitle}{%
 \newpage\pagestyle{empty}
 \begin{tikzpicture}[remember picture,overlay]
  \node[yshift=-\headerheight] at (current page.north west)
   {\begin{tikzpicture}[remember picture, overlay]
     % First draw the grid and then the Logo and avancez logo.
     \draw[black] (0,\headerheight) -- +(\paperwidth,0) -- (\paperwidth,0) -- (0.75\paperwidth,0)
                 .. controls (0.725\paperwidth,0) and (0.725\paperwidth,\headeroffset) ..
                 (0.70\paperwidth,\headeroffset) -- (0,\headeroffset) -- cycle; % (0.5\paperwidth,\headeroffset) to[out=0,in=180] (0.6\paperwidth,0cm);
     %\draw[step=5mm,black] (0,0) grid (\paperwidth,6.4cm);
     %\fill[fill=phd_thesis_blue] (0,0) rectangle (\paperwidth,6.4cm);
     \definecolor{logo_color}{gray}{0}
     \node[anchor=west,xshift=0.1\paperwidth,yshift=0.511\headerheight]{\resizebox{0.63\paperwidth}{!}{\input{figures/CHALMERS.pgf}}}; % 13.2cm
     \node[anchor=west,xshift=0.745\paperwidth,yshift=0.516\headerheight]{\includegraphics[width=0.171\paperwidth]{% 3.6cm
        %figures/Avancez_bg.pdf%
        %figures/Avancez_clean_white.pdf%
        %figures/Avancez_clean_gold.pdf%
        %figures/Avancez_clean_invert.pdf%
        figures/Avancez_clean_black.pdf%
     }};
   \end{tikzpicture}};
 \end{tikzpicture}
 \ \vfill
 \begin{center}
  \@ifundefined{@coverfigure}{}{\includegraphics[width=\textwidth,height=0.4\paperheight,keepaspectratio]{\@coverfigure}}
 \end{center}
 {\huge\noindent\@title\par} % title, 21
 \@ifundefined{@subtitle}{}{{\large\noindent\@subtitle\par}} % subtitle, 16
 {\normalsize\noindent\textit{\@thesis}\par} % 14 thesis in programmename (only bachelor and master's)
 \vskip 5mm
 {\renewcommand{\and}{\\}\Large\noindent \MakeUppercase{\@author}\par}% 18
 \vskip 6mm
 {\small\noindent\@department\\
  \@ifundefined{@division}{}{\textit{\@division}\\}
  \MakeUppercase{\@university}\\
  \@city, \@country\ \@writtenyear\\
  \@ifundefined{@YYYYNN}{}{\@thesis\ \@YYYYNN}\par} % TODO if YYYYNN? Also, thesis in ...

 % Second page, blank
 \newpage
 \mbox{}

 % Third page
 \newpage\thispagestyle{empty}
 \begin{center}
  \@ifundefined{@YYYYNN}{}{{\MakeUppercase{\@thesis}\ \@YYYYNN\par}} % Only if YYYYNN is defined TODO:'thesis in'
  \vspace*{\stretch{1}}%\vskip 40mm
  {\Large \@title\par}
  \vskip 5mm
  \@ifundefined{@subtitle}{}{{\@subtitle\par}}
  \vskip 5mm
  {\@thesis\par} % (if applicable) (if applicable)
  {\renewcommand{\and}{, }\MakeUppercase{\@author}\par}
  \vspace*{\stretch{3}}%\vfill
  {\@department\par}
  \@ifundefined{@division}{}{{\textit{\@division}\par}}
  {\MakeUppercase{\@university}\par}
  \vskip 2mm
  {\@city, \@country\ \@writtenyear\par}
 \end{center}

 % Fourth page
 \newpage
 {\noindent \@title\\
  \@ifundefined{@subtitle}{}{\@subtitle\\}
  \renewcommand{\and}{\\}\MakeUppercase{\@author}
  %\@ifundefined{\@ISBN}{\\ISBN \@ISBN\\}{}
  \par}
 \vskip 10mm
 {\noindent \copyright {\renewcommand{\and}{, }\MakeUppercase{\@author}}, \@writtenyear\par}
 \vspace*{\stretch{1}} %\vskip 20mm
 {\noindent \@ifundefined{@YYYYNN}{@YYYYNN}{}{\@thesis\ \@YYYYNN\\} % TODO: 'thesis in'
  ISSN \@ISSN\\
  \@department\\
  \@ifundefined{@division}{}{\@division\\}
  \@university\\
  \@postcode\ \@city\\
  \@country\\
  \iflanguage{swedish}{Telefon}{Telephone}: \@telephone\par}
 \vspace*{\stretch{3}}
 \@ifundefined{@covercaption}{}{{\noindent \iflanguage{swedish}{Omslag}{Cover}:\\ \@covercaption\par}} % Another translation here
 \vskip 5mm
 {\noindent \@printers
  %\@department
  \\ \@city, \@country\ \@writtenyear\par}

 % Fifth page
 \newpage%\pagestyle{justpage} % TODO page style!
 \setcounter{page}{1}
 \pagenumbering{Roman}
 {\noindent \@title\\
  %\makeatletter\@ifundefined{mysubtitle}{}{\mysubtitle\\}\makeatother
  \@ifundefined{@subtitle}{}{\@subtitle\\}
  \@thesis\\ % (if applicable) thesis in programmename
  {\renewcommand{\and}{\\}\MakeUppercase{\@author}}\\
  \@department\\
  \@ifundefined{@division}{}{\@division\\}
  \@university\par}
 \phantomsection\addcontentsline{toc}{section}{Abstract}\begin{abstract}\@myabstract\end{abstract}
 \noindent \iflanguage{swedish}{Nyckelord}{Keywords}: \@keywords

}

