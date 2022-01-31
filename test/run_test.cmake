set(ARGS_LIST ${ARGS})
separate_arguments(ARGS_LIST)

execute_process(COMMAND ${CMD} ${ARGS_LIST} RESULT_VARIABLE CMD_RESULT COMMAND_ECHO STDOUT)
if(CMD_RESULT)
  message(FATAL_ERROR "Error running ${CMD}: ${CMD_RESULT}")
endif()

execute_process(COMMAND md5sum ${LOG_FILE})
execute_process(COMMAND md5sum ${OUT_FILE})