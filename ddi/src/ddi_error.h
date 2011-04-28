#ifndef DDI_ERROR_H
#define DDI_ERROR_H

/**
   @file
   @brief DDI error macros and operations.
*/

#define DDI_STDERR stderr  /**< @brief DDI standard error */

#define DDI_MAX_PROCESSORS_ERROR    1000  /**< @brief MAX_PROCESSORS exceeded error. */
#define DDI_MAX_NODES_ERROR         1001  /**< @brief MAX_NODES exceeded error. */
#define DDI_MAX_SMP_PROCS_ERROR     1002  /**< @brief MAX_SMP_PROCS exceeded error. */

/**
   @brief Prints the error signal and the optional error message and calls DDI_Abort(int signal). 
   @note If message != NULL, error followed by the error message is printed.
   @note Operating system signals and DDI errors defined in this file have predefined
   error messages and should be called with NULL message.
   @param[in] signal Error signal.
   @param[in] message Optional error message.
 */
void DDI_Error(int signal, char *message);

/**
   @brief Aborts and cleans up DDI runtime.
   @param[in] signal Error signal.
 */
void DDI_Abort(int signal);

#endif /* DDI_ERROR_H */ 
