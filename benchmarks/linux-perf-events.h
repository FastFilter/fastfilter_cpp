// https://github.com/WojciechMula/toys/blob/master/000helpers/linux-perf-events.h
#pragma once
#ifdef __linux__

#include <asm/unistd.h>       // for __NR_perf_event_open
#include <linux/perf_event.h> // for perf event constants
#include <sys/ioctl.h>        // for ioctl
#include <unistd.h>           // for syscall

#include <cerrno>  // for errno
#include <cstring> // for memset
#include <stdexcept>

#include <vector>

template <int TYPE = PERF_TYPE_HARDWARE> class LinuxEvents {
  int fd;
  perf_event_attr attribs;
  int num_events;
  std::vector<uint64_t> temp_result_vec;
  std::vector<uint64_t> ids;

public:
  LinuxEvents(std::vector<int> config_vec) : fd(0) {
    memset(&attribs, 0, sizeof(attribs));
    attribs.type = TYPE;
    attribs.size = sizeof(attribs);
    attribs.disabled = 1;
    attribs.exclude_kernel = 1;
    attribs.exclude_hv = 1;

    attribs.sample_period = 0;
    attribs.read_format = PERF_FORMAT_GROUP | PERF_FORMAT_ID;
    const int pid = 0;  // the current process
    const int cpu = -1; // all CPUs
    const unsigned long flags = 0;

    int group = -1; // no group
    num_events = config_vec.size();
    uint32_t i = 0;
    for (auto config : config_vec) {
      attribs.config = config;
      fd = syscall(__NR_perf_event_open, &attribs, pid, cpu, group, flags);
      if (fd == -1) {
        report_error("perf_event_open");
      }
      ioctl(fd, PERF_EVENT_IOC_ID, &ids[i++]);
      if (group == -1) {
        group = fd;
      }
    }

    temp_result_vec.resize(num_events * 2 + 1);
  }

  ~LinuxEvents() {
      // This solution expects that the fds are increasing numbers without
      // gaps. Certainly this is not always true.
      for (int i = fd; fd - num_events - i != 0; i--) { close(i); }
  }

  inline void start() {
    if (ioctl(fd, PERF_EVENT_IOC_RESET, PERF_IOC_FLAG_GROUP) == -1) {
      report_error("ioctl(PERF_EVENT_IOC_RESET)");
    }

    if (ioctl(fd, PERF_EVENT_IOC_ENABLE, PERF_IOC_FLAG_GROUP) == -1) {
      report_error("ioctl(PERF_EVENT_IOC_ENABLE)");
    }
  }

  inline void end(std::vector<unsigned long long> &results) {
    if (ioctl(fd, PERF_EVENT_IOC_DISABLE, PERF_IOC_FLAG_GROUP) == -1) {
      report_error("ioctl(PERF_EVENT_IOC_DISABLE)");
    }

    if (read(fd, &temp_result_vec[0], temp_result_vec.size() * 8) == -1) {
      report_error("read");
    }
    // our actual results are in slots 1,3,5, ... of this structure
    // we really should be checking our ids obtained earlier to be safe
    for (uint32_t i = 1; i < temp_result_vec.size(); i += 2) {
      results[i / 2] = temp_result_vec[i];
    }
  }

private:
  void report_error(const std::string &context) {
    throw std::runtime_error(context + ": " + std::string(strerror(errno)));
  }
};
#endif
