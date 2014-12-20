#ifndef RUN_CMD_H
#define RUN_CMD_H

#include "logger.hpp"
#include <string>
#include <cstdlib>

/**
 * Runs a system command
 *
 * @param cmd system command
 *
 * @return -1 if command processor is not available, process status
 * otherwise
 */
int run_cmd(const std::string& cmd)
{
   if (!system(NULL)) {
      ERROR("Error: command processor not available!");
      return -1;
   }

   const int status = system(cmd.c_str());

   if (status) {
      VERBOSE_MSG("run_cmd: command \"" << cmd
                  << "\" returned with exit code " << status);
   }

   return status;
}

#endif
