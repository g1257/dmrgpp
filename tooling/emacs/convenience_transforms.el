(defun my/transform-dmrg-tests-to-awk-check (start end)
  "Transform a region of add_test lines for dmrg into add_awk_test lines.
For each line like:
  add_test(NAME inputN COMMAND ./dmrg -f ${PATH_TO_INPUTS}/inputN.ain) # energy VALUE
Insert a new line after it:
  add_lowest_eigenvalue_check_test(outputN VALUE runForinputN.cout)
Has minor bug that when 3 or more lines in region the last does not have
its check added, just select it by itself to work around for now"
  (interactive "r")
  (save-excursion
    (goto-char start)
    (let ((re (rx
               "add_test(NAME"
               (zero-or-more space)
               (group-n 1 (one-or-more (not (any whitespace))))
               (zero-or-more space)
               "COMMAND"
               (zero-or-more space)
               "./dmrg"
               (zero-or-more space)
               "-f"
               (zero-or-more space)
               "${PATH_TO_INPUTS}/"
               (group-n 2 (one-or-more (not (or (any ?/) (any whitespace)))))
               ".ain)"
               (zero-or-more space)
               "# energy"
               (zero-or-more space)
               (group-n 3 (optional "-")
                        (one-or-more digit)
                        "."
                        (one-or-more digit)))))
      (while (< (point) (+ 1 end))
        (let ((line (buffer-substring-no-properties
                     (point)
                     (progn (end-of-line) (point)))))
          (message "re = %s" re)
          (message "line = %S" line)
          (message "string-match result: %s" (string-match re line))

          (when (string-match re line)
            (let* ((name (match-string 1 line))
                   (basename (match-string 2 line))
                   (expected (match-string 3 line))
                   (new-name (replace-regexp-in-string "input" "output" name))
                   (output-file (concat "runFor" basename ".cout")))
              (insert "\nadd_lowest_eigenvalue_check_test(" new-name " " expected " " output-file ")\n")))
          (forward-line 1))))))
