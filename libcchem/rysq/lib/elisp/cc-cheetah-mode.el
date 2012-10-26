(require 'cc-mode)

;; These are only required at compile time to get the sources for the
;; language constants.  (The cc-fonts require and the font-lock
;; related constants could additionally be put inside an
;; (eval-after-load "font-lock" ...) but then some trickery is
;; necessary to get them compiled.)
(eval-when-compile
  (require 'cc-langs)
  (require 'cc-fonts))

(eval-and-compile
  ;; Make our mode known to the language constant system.  Use Java
  ;; mode as the fallback for the constants we don't change here.
  ;; This needs to be done also at compile time since the language
  ;; constants are evaluated then.
  (c-add-language 'c++-cheetah-mode 'c++-mode))
(c-lang-defconst c-opt-cpp-start c++-cheetah "\\s *\\(#\\|%\\)\\s *\\([[:alnum:]]+\\)" )

(defcustom c++-cheetah-font-lock-extra-types nil
  "*List of extra types (aside from the type keywords) to recognize in C++-Cheetah mode.
Each list item should be a regexp matching a single identifier.")

(defconst c++-cheetah-font-lock-keywords-1 (c-lang-const c-matchers-1 c++-cheetah)
  "Minimal highlighting for C++-Cheetah mode.")

(defconst c++-cheetah-font-lock-keywords-2 (c-lang-const c-matchers-2 c++-cheetah)
  "Fast normal highlighting for C++-Cheetah mode.")

(defconst c++-cheetah-font-lock-keywords-3 (c-lang-const c-matchers-3 c++-cheetah)
  "Accurate normal highlighting for C++-Cheetah mode.")

(defvar c++-cheetah-font-lock-keywords c++-cheetah-font-lock-keywords-3
  "Default expressions to highlight in blah mode.")

(font-lock-add-keywords 'c++-cheetah-mode
  '(("%\\(end\\ [a-z]+\\|[a-z]+\\)" 0 font-lock-preprocessor-face t)
    ("\\$(\\w+)" 0 font-lock-preprocessor-face t)
    ("%.+\\ *\\(#.*\\)" 1 font-lock-comment-face t)))

(defvar c++-cheetah-mode-syntax-table nil
  "Syntax table used in c++-cheetah-mode buffers.")
(or c++-cheetah-mode-syntax-table
    (setq c++-cheetah-mode-syntax-table
	  (funcall (c-lang-const c-make-mode-syntax-table c++-cheetah))))

;; (defvar c++-cheetah-mode-abbrev-table nil
;;   "Abbreviation table used in c++-cheetah-mode buffers.")
;; (c-define-abbrev-table 'c++-cheetah-mode-abbrev-table
;;   ;; Keywords that if they occur first on a line might alter the
;;   ;; syntactic context, and which therefore should trig reindentation
;;   ;; when they are completed.
;;   '(("else" "else" c-electric-continued-statement 0)
;;     ("while" "while" c-electric-continued-statement 0)
;;     ("catch" "catch" c-electric-continued-statement 0)
;;     ("finally" "finally" c-electric-continued-statement 0)))



(defvar c++-cheetah-mode-map (let ((map (c-make-inherited-keymap)))
		      ;; Add bindings which are only useful for C++-Cheetah
		      map)
  "Keymap used in c++-cheetah-mode buffers.")

;; (easy-menu-define c++-cheetah-menu c++-cheetah-mode-map "C++-Cheetah Mode Commands"
;; 		  ;; Can use `c++-cheetah' as the language for `c-mode-menu'
;; 		  ;; since its definition covers any language.  In
;; 		  ;; this case the language is used to adapt to the
;; 		  ;; nonexistence of a cpp pass and thus removing some
;; 		  ;; irrelevant menu alternatives.
;; 		  (cons "C++-Cheetah" (c-lang-const c-mode-menu c++-cheetah)))

;;;###autoload
(add-to-list 'auto-mode-alist '("\\.\\(cpp\\|c\\|h\\)\\.tmpl\\'" . c++-cheetah-mode))

;;;###autoload
(defun c++-cheetah-mode ()
  "Major mode for editing C++-Cheetah (pronounced \"big nose\") code.
This is a simple example of a separate mode derived from CC Mode to
support a language with syntax similar to C/C++/ObjC/Java/IDL/Pike.

The hook `c-mode-common-hook' is run with no args at mode
initialization, then `c++-cheetah-mode-hook'.

Key bindings:
\\{c++-cheetah-mode-map}"
  (interactive)
  (kill-all-local-variables)
  (c-initialize-cc-mode t)
  (set-syntax-table c++-cheetah-mode-syntax-table)
  (setq major-mode 'c++-cheetah-mode
	mode-name "C++-Cheetah"
;;	local-abbrev-table c++-cheetah-mode-abbrev-table
;;	abbrev-mode t
	)
  (use-local-map c-mode-map)
  ;; `c-init-language-vars' is a macro that is expanded at compile
  ;; time to a large `setq' with all the language variables and their
  ;; customized values for our language.
  (c-init-language-vars c++-cheetah-mode)
  ;; `c-common-init' initializes most of the components of a CC Mode
  ;; buffer, including setup of the mode menu, font-lock, etc.
  ;; There's also a lower level routine `c-basic-common-init' that
  ;; only makes the necessary initialization to get the syntactic
  ;; analysis and similar things working.
  (c-common-init 'c++-cheetah-mode)
;;  (easy-menu-add c++-cheetah-menu)
  (run-hooks 'c-mode-common-hook)
  (run-hooks 'c++-cheetah-mode-hook)
  (c-update-modeline))

(provide 'c++-cheetah-mode)

;;; derived-mode-ex.el ends here
