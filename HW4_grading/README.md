# Homework 4 Grading Script
We will use this script to grade your program. **Make sure your program can be executed by this script.**

## Preparing
* Step1:  
    Enter the `HW4_grading` directory and create a new directory named with your student ID in the `student` directory.
    ```sh
    $ cd HW4_grading/
    $ mkdir student/${your_student_id}
    ```
    For example,
    ```sh
    $ cd HW4_grading/
    $ mkdir student/112062500
    ```

* Step2:  
    Put your compressed file in the directory which you just created.  
    The correct path should be:
    ```
    HW4_grading/student/${your_student_id}/CS6135_HW4_${your_student_id}.tar.gz
    ```
    For example,
    ```
    HW4_grading/student/112062500/CS6135_HW4_112062500.tar.gz
    ```

### Notice:
**Please make sure not to put your original directory here**, as it will remove all directories before unzipping the compressed file.

## Grading
* Step1:  
    Navigate to the `HW4_grading` directory and run `HW4_grading.sh`.
    ```sh
    $ cd HW4_grading/
    $ bash HW4_grading.sh
    ```

* Step2:  
    Check your output. Ensure the status of each checking item is **yes**.
    * If the status of a testcase is **success**, it means your program finished in time, and the output result is legal.
        ```
        grading on 112062500:
         checking item          | status
        ------------------------|--------
         correct tar.gz         | yes
         correct file structure | yes
         have README            | yes
         have Makefile          | yes
         correct make clean     | yes
         correct make           | yes

          testcase | wirelength |    runtime | status
        -----------|------------|------------|--------
           public1 |  324193780 |       0.18 | success
           public2 |   29503287 |       0.44 | success
           public3 | 2644255425 |       0.79 | success
        ```

    * If the status of a testcase is not **success**, it means your program failed in this testcase.
        ```
        grading on 112062500:
         checking item          | status
        ------------------------|--------
         correct tar.gz         | yes
         correct file structure | yes
         have README            | yes
         have Makefile          | yes
         correct make clean     | yes
         correct make           | yes
        
          testcase | wirelength |    runtime | status
        -----------|------------|------------|--------
           public1 |        N/A |        N/A | There is an error in the output results of public1 (...).
           public2 |        N/A |        N/A | There is an error in the output results of public2 (...).
           public3 |        N/A |        N/A | There is an error in the output results of public3 (...).
        ```
