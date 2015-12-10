% --- identification
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{standard-doc}[2013/03/10 note class]

% --- using pkgs
\RequirePackage{fancyhdr}
\RequirePackage{lastpage}


% --- define options
\newif\if@color\@colorfalse
\newif\if@pdf\@pdffalse
\DeclareOption{color}{\@colortrue}
\DeclareOption{pdftex}{\@pdftrue}
\DeclareOption{twocolumn}{\OptionNotUsed}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax

\if@pdf
\RequirePackage[pdftex]{graphicx}
\else
\RequirePackage{graphicx}
\fi

\LoadClass[a4paper]{article}

\topmargin -15mm
\headsep 10mm
\headheight 23mm

%%\footheight 0mm

\textheight230mm 
\textwidth170mm

\hoffset 0mm
\voffset 0mm
%  modif centrage Coulon (papier A4, pas Letter) :
\oddsidemargin  -5mm
\evensidemargin -5mm

\leftmargin 20mm
\rightmargin 20mm


% --- load class


% color is switched on, use the color package
\if@color
\RequirePackage{color}
\fi%color




\def\issue#1{\gdef\@issue{#1}}
\def\@issue{\@latex@error{No \noexpand\issue given}\@ehc}
\def\shorttitle#1{\gdef\@shorttitle{#1}}
\def\@shorttitle{\@latex@error{No \noexpand\shorttitle given}\@ehc}

\renewcommand\maketitle{
\begin{titlepage}\begin{center}

%\begin{minipage}[t]{.3\textwidth}
%{\includegraphics[height=2.0cm]{bhamcrest.png}} 
%\end{minipage}
%\begin{minipage}[t]{.3\textwidth}
%{\includegraphics[height=1.5cm]{bham_word_marque.png}} 
%\end{minipage}

\vskip2cm
   {\LARGE \bf \@title \par}%
   \vskip 3em%
   {\large
     \lineskip .75em%
     \begin{tabular}[t]{c}%
       \@author
     \end{tabular}\par}%
		 \ifx\@ligodcc\@empty
		 \else
		 \vskip 1.5em%
		 \fi%ligodcc
		 \vskip 1.5em%
   {\normalsize \textit{Issue:}~ \@issue \par}%
   {\normalsize \textit{Date:}~ \@date \par}%
   \vfil\null
\end{center}
\end{titlepage}
}


\pagestyle{fancy}
%\lhead{\setlength{\unitlength}{1mm}
%{\includegraphics[height=10mm]{bhamlogo.pdf}}}

\chead{\@shorttitle}
\rhead{}

\lfoot{} 
\cfoot{{\footnotesize \textit{page} \thepage\ of \pageref{LastPage}}} 
\rfoot{}




