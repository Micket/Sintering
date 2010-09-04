%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Eccomas class for LaTeX files                 v1.1 Feb-23-2000 %% <-- Outdated
%%                                                                %%
%% The regular article class commands still exist and may be used.%%
%% 10pt, 11pt and twocolumn options are disabled.                 %%
%%                                                                %%
%% \maketitle make the first page automatically                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\NeedsTeXFormat{LaTeX2e}
%\ProvidesClass{membranes}[2007/09/17-19 MEMBRANES paper class]
\ProvidesClass{nordic}[2008/10/16-17 NORDIC paper class]

\DeclareOption{10pt}{\OptionNotUsed}
\DeclareOption{11pt}{\OptionNotUsed}
\DeclareOption{twocolumn}{\OptionNotUsed}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}

\ProcessOptions\relax

\LoadClass[11pt]{article}

\RequirePackage{calc}
\RequirePackage{indentfirst}
% I think biblatex is still to new for most users, and natbib is well established.
\RequirePackage[super,comma]{natbib} % Turns normal \cite into the required format.

%% Title
\renewcommand{\maketitle}{%
 \newpage
 \thispagestyle{plain}
 \begin{minipage}[t][1cm][t]{15.5cm}
  \begin{flushright}
   \scriptsize
   \vspace{-50pt}
   \@header
  \end{flushright}
 \end{minipage}
 \begin{center}%
  \vspace{-18pt}%
  {%3
   \fontsize{14}{16}\selectfont
   \bfseries
   \MakeUppercase{\@title}
   \par
  }%3
  %\vspace{18pt}%
  \vspace{13pt}%
  {%4
   \fontsize{12}{14}\selectfont
   \bfseries\
   \MakeUppercase{\@author}
   \par
  }%4
  %\vspace{14pt}%
  \vspace{10pt}%
  {%5
   \fontsize{11}{12}\selectfont
   %\def\and{\\\vspace{14pt}}
   \def\and{\\\vspace{13pt}}
   \normalfont\@address
   \par
  }%5
 \end{center}
 %\vspace{18pt}%
 \vspace{11pt}%
 \noindent\textbf{Key words:\space}\@keywords
 \par
 %\vspace{18pt}%
 \vspace{12pt}%
 \noindent\textbf{Summary.\space}{\normalfont\@abstract}
 %\noindent\textbf{Abstract.\space}{\textit\@abstract}
 \par
 \vspace{7pt}% Removes the annoying \vskip-0.4cm in the document
 %\vfill
 %\newpage
}% End of \@maketitle

\def\author{\@ifnextchar[{\@getshortauthor}{\@getauthor}}
\def\@getauthor#1{\gdef\@author{#1}\gdef\@shortauthor{#1}}
\def\@getshortauthor[#1]#2{\gdef\@shortauthor{#1}\gdef\@author{#2}}

\newcommand{\header}[1]{\gdef\@header{#1}}

\newcommand{\address}[1]{\gdef\@address{#1}}

\newcommand{\keywords}[1]{\gdef\@keywords{#1}}

\renewcommand{\abstract}[1]{\gdef\@abstract{#1}}

\newcommand{\heading}[1]{\gdef\@heading{#1}}

\def\ps@myheadings{%
    \def\@oddfoot{\normalfont\hfil\thepage\hfil}
    \let\@evenfoot\@oddfoot
    \def\@evenhead{%
      %\parbox[t][11.5pt][t]{\textwidth}{%
      \parbox[t][16.5pt][t]{\textwidth}{%
        \centering\normalfont\fontsize{10}{10}\selectfont\@heading\\
        \vspace{-.2cm}%
        \rule{\textwidth}{1pt}%
     }
    }%
    \let\@oddhead\@evenhead
    \let\@mkboth\@gobbletwo
    \let\sectionmark\@gobble
    \let\subsectionmark\@gobble
}

\def\sectionLevel{1}
\def\sectionIndent{0pt}
\def\sectionSpaceBefore{12pt}
\def\sectionSpaceAfter{6pt}
\def\sectionStyle{\normalsize\bfseries\MakeUppercase}

\renewcommand{\section}{%
  \@startsection{section}{\sectionLevel}{\sectionIndent}{\sectionSpaceBefore}%
                {\sectionSpaceAfter}{\sectionStyle}}%

\def\sectionLevel{2}
\def\sectionIndent{0pt}
\def\sectionSpaceBefore{12pt}
\def\sectionSpaceAfter{6pt}
\def\sectionStyle{\normalsize\bfseries}

\renewcommand{\subsection}{%
  \@startsection{subsection}{\sectionLevel}{\sectionIndent}{\sectionSpaceBefore}%
                {\sectionSpaceAfter}{\sectionStyle}}%

\renewcommand{\@makecaption}[2]{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{\normalsize\fontsize{10}{12}\selectfont #1: #2}%
  \ifdim \wd\@tempboxa >\hsize
    \normalsize\fontsize{10}{12}\selectfont #1: #2\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}

\renewenvironment{thebibliography}[1]
  {%\newpage
   \section*{\MakeUppercase\refname\@mkboth{\MakeUppercase\refname}{\MakeUppercase\refname}}%
   \list{\@biblabel{\@arabic\c@enumiv}}%
        {\settowidth\labelwidth{\@biblabel{#1}}%
         \leftmargin\labelwidth
         \advance\leftmargin\labelsep
         \@openbib@code
         \usecounter{enumiv}%
         \let\p@enumiv\@empty
         \renewcommand\theenumiv{\@arabic\c@enumiv}}%
   \sloppy
   \clubpenalty4000
   \@clubpenalty \clubpenalty
   \widowpenalty4000%
   \sfcode`\.\@m}
  {\def\@noitemerr{\@latex@warning{Empty `thebibliography' environment}}%
  \endlist}
\setlength{\oddsidemargin}{0mm}
\setlength{\textheight}{210mm}
\setlength{\textwidth}{160mm}
\setlength{\intextsep}{12pt} % This is what you meant to set to 12pt, not textfloatsep
\setlength{\abovedisplayskip}{6pt}
\setlength{\belowdisplayskip}{6pt}
\setlength{\parindent}{5mm}
\setlength{\abovecaptionskip}{6pt}
\setlength{\belowcaptionskip}{0pt} % Actual skip = \textfloatsep + \belowcaptionskip
\setlength{\itemsep}{0pt}

\AtBeginDocument{%
  \pagestyle{myheadings}
  %\maketitle % In order to be able to use inputencodings and urls, this can't be done here.
  %\let\maketitle\relax
}

\endinput
