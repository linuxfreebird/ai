
(dribble "output2.text")
(print "hello world")

(defun double (x) (* x 2))
(print (double 2))
(DEFUN factorial (N)
  "Compute the factorial of N."
  (if (= N 1)
      1
    (* N (factorial (- N 1)))))
(print (factorial 3))

(defun fibonacci (N)
  "Compute the N'th Fibonacci number."
  (if (or (zerop N) (= N 1))
      1
    (let
        (
         (F1 (fibonacci (- N 1)))
         (F2 (fibonacci (- N 2)))
         )
      (+ F1 F2)
      )
    )
  )

(print (fibonacci 5))
