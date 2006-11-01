
#ifndef _libint2_src_bin_libint_task_h_
#define _libint2_src_bin_libint_task_h_

#include <string>
#include <vector>
#include <default_params.h>

namespace libint2 {

  /**
     A key idea introduced here is that of "task". By task I mean a
     type of a computation that Libint performs. For example, computing ERIs is one task,
     computing sets of integrals needed for the R12 theory is another task. The reason for speaking of separate tasks
     is that the evaluator type must be specific to each task.
     For example, all the intermediates needed in R12 computation are not necessary when computing ERIs.
     If one evaluator type were to cover all tasks, it would be huge and performance would be likely hurt.
     Thus we need to produce task-specific evaluator and associated functions (constructor, destructor,
     memory query).
  */
  class LibraryTask {
  public:
    LibraryTask(const std::string& l, const SafePtr<TaskParameters>& p) : label_(l), params_(p) {
    }
    ~LibraryTask() {
    }
    const std::string& label() const { return label_; }
    const SafePtr<TaskParameters>& params() const { return params_; }

  private:
    std::string label_;
    SafePtr<TaskParameters> params_;
  };

  /// Manages tasks. This is a Singleton.
  class LibraryTaskManager {
    typedef std::vector<LibraryTask> Tasks;
    typedef Tasks::const_iterator tasks_citer;

  public:
    typedef Tasks::const_iterator TasksCIter;
    typedef Tasks::iterator TasksIter;

    static LibraryTaskManager& Instance();
    ~LibraryTaskManager() {}

    /// Number of tasks
    unsigned int ntask() const { return tasks_.size(); }
    /// returns iterator to the first task
    TasksCIter first() const { return tasks_.begin(); }
    /// returns iterator to past the last task
    TasksCIter last() const { return tasks_.end(); }
    /// i-th tasks
    LibraryTask& task(unsigned int i) { return tasks_.at(i); }
    /// Adds a new task
    void add(const std::string& task_label);
    /// Makes this task current (must have been added already)
    void current(const std::string& task_label);
    /// Returns the current task
    LibraryTask& current();

  private:

    LibraryTaskManager();
    Tasks tasks_;
    int current_;

    static LibraryTaskManager LTM_obj_;
  };

};

#endif // header guard
