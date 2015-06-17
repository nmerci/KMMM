%%
%% LaTeX documentclass for the journal
%% Теорія ймовірностей та математична статистика
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%<идентификация>%
\NeedsTeXFormat{LaTeX2e}
	[1994/12/01]
\ProvidesClass{tbimc}
	[2004/03/31 v1.0471 TBIMC Author Class]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%<начальный командный код>%
% Define commands for project definition file

% Process page number as a global counter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%<декларация опций>%
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{amsart}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%<выполнение опций>%
\ProcessOptions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%<загрузка пакетов>%
\RequirePackage[T2A]{fontenc}
\RequirePackage[cp1251]{inputenc}
\RequirePackage[english,russian,ukrainian]{babel}

\LoadClass[reqno]{amsart}[1999/11/24]

\RequirePackage{amssymb}

\textwidth=32pc
\textheight=52.5pc
\calclayout

%% Change text "and" to \andname
\addto\extrasukrainian{\renewcommand{\andify}{%
  \nxandlist{\unskip, }{\unskip{} \andname~}{\unskip{} \andname~}}}
\addto\extrasrussian{\renewcommand{\andify}{%
  \nxandlist{\unskip, }{\unskip{} \andname~}{\unskip{} \andname~}}}
\addto\extrasenglish{\renewcommand{\andify}{%
  \nxandlist{\unskip, }{\unskip{} \andname~}{\unskip, \andname~}}}

%% Some translations for captions
\addto\captionsukrainian{%
	\def\UDCname{УДК}%
	\def\keywordsname{Ключові слова і фрази}%
	\def\datename{Надійшла}%
	\def\curraddrname{Поточна адреса}%
	\def\emailname{Адреса електронної пошти}%
	\def\andname{та}%
}
\addto\captionsrussian{%
	\def\UDCname{УДК}%
	\def\keywordsname{Ключевые слова и фразы}%
	\def\datename{Поступила}%
	\def\curraddrname{Текущий адрес}%
	\def\emailname{Адрес электронной почты}%
	\def\andname{и}%
}
\addto\captionsenglish{%
	\def\UDCname{UDC}%
	\def\keywordsname{Key words and phrases}%
	\def\datename{Received}%
	\def\curraddrname{Current address}%
	\def\emailname{E-mail address}%
	\def\andname{and}%
}

%% Default language
\selectlanguage{ukrainian}

%% Set lower and upper case for Ukrainian letters in cp1251 encoding
\uccode`\ґ=`\Ґ	\lccode`\Ґ=`\ґ
\uccode`\є=`\Є	\lccode`\Є=`\є
\uccode`\і=`\І	\lccode`\І=`\і
\uccode`\ї=`\Ї	\lccode`\Ї=`\ї
%%

%% Set margins
\setlength{\oddsidemargin}{0pt}
\setlength{\evensidemargin}{0pt}
\setlength{\topmargin}{-2pt}
\addtolength{\textheight}{10pt}

%% Sider

%%
\endinput