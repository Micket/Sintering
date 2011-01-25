% Thesis class template for Chalmers University of Technology.
% This is an inofficial class made by Mikael Öhman <mikael.ohman@chalmers.se>
% This document is released under public domain
% Some words are translated to swedish if babel langauge is set, otherwise default to English

\NeedsTeXFormat{LaTeX2e}

\ProvidesClass{chalmers-thesis}[2010/09/09 CHALMERS Thesis class]

%\LoadClass[12pt]{scrartcl}
\LoadClass[a4paper,12pt]{article}

\newif\if@bachelor
\newif\if@masters
\newif\if@licenciate
\newif\if@doctorate

\@bachelorfalse
\@mastersfalse
\@licenciatefalse
\@doctoratefalse

%% Options of actamont class
\DeclareOption{bachelor}{\global\@bachelortrue
  \global\@mastersfalse\global\@licenciatefalse\global\@doctoratefalse}
\DeclareOption{masters}{\global\@masterstrue
  \global\@bachelorfalse\global\@licenciatefalse\global\@doctoratefalse}
\DeclareOption{licenciate}{\global\@licenciatetrue
  \global\@bachelorfalse\global\@mastersfalse\global\@doctoratefalse}
\DeclareOption{doctorate}{\global\@doctoratetrue
  \global\@bachelorfalse\global\@masterstrue\global\@licenciatefalse}

%\DeclareOption{bachelor}{\def\@bachelortrue{}}
%\DeclareOption{masters}{\def\@masterstrue{}}
%\DeclareOption{licenciate}{\def\@licenciatetrue{}}
%\DeclareOption{doctorate}{\def\@doctoratetrue{}}

%\DeclareOption*{\PassOptionsToClass{\CurrentOption}{scrartcl}}
%\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}

% execute default options
\ExecuteOptions{masters}
\ProcessOptions

\RequirePackage{tikz}
\RequirePackage{fancyhdr}
\RequirePackage{hyperref}
\RequirePackage{graphicx}
\RequirePackage[swedish,english]{babel}

\usepackage[top=0.1\paperheight, bottom=0.1\paperheight, left=0.1\paperwidth, right=0.1\paperwidth]{geometry}

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
%\renewcommand{\title}[1]{\gdef\@title{#1}}
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
\country{\iflanguage{swedish}{Sverige}{Sweden}}
\university{\iflanguage{swedish}{Chalmers Tekniska Högskola}{Chalmers University of Technology}}
\university{Chalmers University of Technology} % expandafter just won't work, MakeUppercase is annoying.
\printers{Chalmers Reproservice}
\telephone{+ 46 (0)31-772 1000}
\if@bachelor
  \ISSN{1654-4676} % Bachelor
  \def\titlepage@edge{draw}
  \def\titlepage@chalmerscolor{black}
  \def\titlepage@bg{white}
  \def\titlepage@logotype{Avancez_clean_invert}
\fi
\if@masters
  \ISSN{1652-8557} % Masters
  \def\titlepage@edge{draw}
  \def\titlepage@chalmerscolor{black}
  \def\titlepage@bg{step=5mm}
  \def\titlepage@logotype{Avancez_clean_white}
\fi
\if@licenciate
  \ISSN{1652-8565} % Lic.
  \def\titlepage@edge{}
  \def\titlepage@chalmerscolor{white}
  \def\titlepage@bg{fill=chalmers_grey}
  \def\titlepage@logotype{Avancez_clean_black}
\fi
\if@doctorate
  \ISSN{0346-718X} % PhD
  \def\titlepage@edge{}
  \def\titlepage@chalmerscolor{white}
  \def\titlepage@bg{fill=chalmers_blue}
  \def\titlepage@logotype{Avancez_clean_gold}
\fi

% Introducing some new lengths for the front page
\newlength{\headerheight}
\newlength{\headeroffset}
\setlength{\headerheight}{0.215\paperheight}%~6.4cm
\setlength{\headeroffset}{0.125\headerheight}%~0.8cm

% Document metadata
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
     \clip[\titlepage@edge] (0,\headerheight) -- +(\paperwidth,0) -- (\paperwidth,0) -- (0.75\paperwidth,0)
                 .. controls (0.725\paperwidth,0) and (0.725\paperwidth,\headeroffset) ..
                 (0.70\paperwidth,\headeroffset) -- (0,\headeroffset) -- cycle; % (0.5\paperwidth,\headeroffset) to[out=0,in=180] (0.6\paperwidth,0cm);
     %\draw[\titlepage@bg](0,0) rectangle (\paperwidth,6.4cm);
     \draw[step=5mm,black] (0,0) grid (\paperwidth,6.4cm);
     %\fill[fill=chalmers_grey] (0,0) rectangle (\paperwidth,6.4cm);
     \node[anchor=west,xshift=0.1\paperwidth,yshift=0.511\headerheight]{\chalmerslogo{\titlepage@chalmerscolor}{0.63\paperwidth}}; % 13.2cm for A4 paper
     \node[anchor=west,xshift=0.745\paperwidth,yshift=0.516\headerheight]{\includegraphics[width=0.171\paperwidth]{% 3.6cm
        figures/\titlepage@logotype
     }};
   \end{tikzpicture}};
 \end{tikzpicture}
 \ \vfill
 \begin{center}
  \@ifundefined{@coverfigure}{}{\@coverfigure}
 \end{center}
 {\huge\noindent\@title\par} % title, 21
 \@ifundefined{@subtitle}{}{{\large\noindent\@subtitle\par}} % subtitle, 16
 {\normalsize\noindent\textit{\@thesis}\par} % 14 thesis in programmename (only bachelor and master's)
 \vskip 1em
 {\renewcommand{\and}{\\}\Large\noindent \MakeUppercase{\@author}\par}% 18
 \vskip 1em
 {\small\noindent\@department\\
  \@ifundefined{@division}{}{\textit{\@division}\\}
  \MakeUppercase{\@university}\\
  \@city, \@country\ \@writtenyear\\
  \@ifundefined{@YYYYNN}{}{\@thesis\ \@YYYYNN}\par} % TODO if YYYYNN? Also, thesis in ...

 % Second page, blank
 \newpage
 \mbox{}

 % Third page, another title page
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

 % Fourth page, printing information
 \newpage
 {\noindent \@title\\
  \@ifundefined{@subtitle}{}{\@subtitle\\}
  \renewcommand{\and}{\\}\MakeUppercase{\@author}
  %\@ifundefined{\@ISBN}{\\ISBN \@ISBN\\}{}
  \par}
 \vskip 10mm
 {\noindent \copyright\ {\renewcommand{\and}{, }\MakeUppercase{\@author}}, \@writtenyear\par}
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

 % Fifth page, first abstract
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

 \newpage
 \pagestyle{fancy}
}

% CHALMERS logo below. Note! It is trademarked!
\newcommand{\chalmerslogo}[2]{
\resizebox{#2}{!}{
\begin{tikzpicture}[fill=#1] % This is 11.994 in width, 1.341 in height
% C
\fill (1.481,0.488)
 .. controls (1.385,0.138) and (1.046,0.027) .. (0.742,0.027)
 .. controls (0.304,0.027) and (0,0.256) .. (0,0.687)
 .. controls (0,1.152) and (0.37,1.34) .. (0.76,1.34)
 .. controls (0.901,1.34) and (1.273,1.307) .. (1.437,1.011)
 -- (1.128,0.914)
 .. controls (1.095,0.966) and (0.997,1.102) .. (0.758,1.102)
 .. controls (0.478,1.102) and (0.321,0.919) .. (0.321,0.683)
 .. controls (0.321,0.414) and (0.535,0.271) .. (0.769,0.271)
 .. controls (1.036,0.271) and (1.116,0.411) .. (1.156,0.486)
 -- cycle;

% H
\fill (1.616,1.302) -- (1.925,1.302) -- (1.925,0.827) -- (2.656,0.827) -- (2.656,1.302) -- (2.965,1.302)
 -- (2.965,0.069) -- (2.656,0.069) -- (2.656,0.598) -- (1.925,0.598) -- (1.925,0.069) -- (1.616,0.069)
 -- cycle;

%A
\fill (3.707,1.302) -- (4.041,1.302) -- (4.637,0.069) -- (4.291,0.069) -- (4.168,0.338) -- (3.568,0.338)
 -- (3.449,0.069) -- (3.112,0.069) -- (3.707,1.302) -- cycle
 (3.903,0.923)
 .. controls (3.889,0.955) and (3.876,0.985) .. (3.871,1.021)
 .. controls (3.866,0.986) and (3.851,0.954) .. (3.837,0.922)
 -- (3.668,0.556) -- (4.066,0.556) -- cycle;

% L
\fill (4.796,1.302) -- (5.105,1.302) -- (5.105,0.296) -- (5.728,0.296) -- (5.728,0.069) -- (4.796,0.069) -- cycle;

% M
\fill (5.914,1.302) -- (6.359,1.302) -- (6.687,0.532)
 .. controls (6.721,0.451) and (6.73,0.418) .. (6.743,0.364)
 .. controls (6.756,0.418) and (6.767,0.451) .. (6.8,0.532)
 -- (7.122,1.302) -- (7.564,1.302) -- (7.564,0.069) -- (7.272,0.069) -- (7.272,0.937)
 .. controls (7.272,0.984) and (7.279,1.01) .. (7.285,1.047)
 .. controls (7.272,1.01) and (7.263,0.984) .. (7.243,0.937)
 -- (6.881,0.069) -- (6.6,0.069) -- (6.231,0.937)
 .. controls (6.211,0.984) and (6.205,1.01) .. (6.195,1.047)
 .. controls (6.198,1.01) and (6.205,0.984) .. (6.205,0.937)
 -- (6.205,0.069) -- (5.914,0.069) -- cycle;

% E
\fill (7.85,1.302) -- (8.999,1.302) -- (8.999,1.080) -- (8.159,1.080) -- (8.159,0.802) -- (8.899,0.802)
 -- (8.899,0.584) -- (8.159,0.584) -- (8.159,0.289) -- (9.034,0.289) -- (9.034,0.069) -- (7.850,0.069) -- cycle;

% R
\fill (9.266,1.301) -- (10.128,1.301)
 .. controls (10.219,1.301) and (10.566,1.290) .. (10.566,0.933)
 .. controls (10.566,0.720) and (10.439,0.612) .. (10.246,0.583)
 -- (10.620,0.068) -- (10.265,0.068) -- ( 9.927,0.564) -- ( 9.575,0.564) -- ( 9.575,0.068)
 -- ( 9.266,0.068) -- ( 9.266,1.301) -- cycle
 (9.576,0.782) -- (10.068,0.782)
 .. controls (10.134,0.782) and (10.253,0.791) .. (10.253,0.929)
 .. controls (10.253,1.064) and (10.164,1.085) .. (10.054,1.085)
 -- (9.576,1.085) -- cycle;

% S
\fill (10.956,0.445)
 .. controls (11.028,0.297) and (11.205,0.255) .. (11.373,0.258)
 .. controls (11.448,0.259) and (11.684,0.276) .. (11.681,0.426)
 .. controls (11.680,0.515) and (11.593,0.534) .. (11.527,0.545)
 .. controls (11.476,0.554) and (11.222,0.589) .. (11.164,0.599)
 .. controls (11.051,0.617) and (10.758,0.667) .. (10.758,0.946)
 .. controls (10.758,1.249) and (11.080,1.340) .. (11.326,1.340)
 .. controls (11.612,1.340) and (11.831,1.285) .. (11.981,1.056)
 -- (11.689,0.978)
 .. controls (11.620,1.103) and (11.462,1.121) .. (11.330,1.120)
 .. controls (11.194,1.120) and (11.055,1.089) .. (11.057,0.972)
 .. controls (11.058,0.892) and (11.130,0.872) .. (11.223,0.853)
 .. controls (11.274,0.843) and (11.526,0.805) .. (11.584,0.795)
 .. controls (11.727,0.770) and (11.971,0.723) .. (11.976,0.452)
 .. controls (11.983,0.056) and (11.539,0.033) .. (11.358,0.031)
 .. controls (11.156,0.028) and (10.840,0.057) .. (10.676,0.366)
 -- cycle;
\end{tikzpicture}
}
}