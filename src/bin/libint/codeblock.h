
#include <string>
#include <entity.h>

#ifndef _libint2_src_bin_libint_codeblock_h_
#define _libint2_src_bin_libint_codeblock_h_

using namespace std;

namespace libint2 {
  
  class CodeContext;
  
  class CodeBlock {
    public:
    CodeBlock(const SafePtr<CodeContext>& context) :
      context_(context) {}
    virtual ~CodeBlock() {}
    
    SafePtr<CodeContext> context() const { return context_; }
    
    /// Opens a code block
    virtual std::string open() =0;
    /// Close a code block
    virtual std::string close() =0;
    
    private:
    SafePtr<CodeContext> context_;
  };
  
  class ForLoop : public CodeBlock {
    public:
    ForLoop(const SafePtr<CodeContext>& context, std::string& varname,
            const SafePtr<Entity>& less_than, const SafePtr<Entity>& start_at);
    virtual ~ForLoop();
    
    /// Implementation of CodeBlock::open()
    std::string open();
    /// Implementation of CodeBlock::close()
    std::string close();
    
    private:
    std::string varname_;
    SafePtr<Entity> less_than_;
    SafePtr<Entity> start_at_;
    
    // checks less_than_ and start_at_ and initializes
    // lt_expr_, sa_expr_, and dummy_loop_
    void init_();
    // string representations of less_than_ and start_at_
    std::string lt_expr_;
    std::string sa_expr_;
    // dummy_loop_ is true if the loop is guaranteed to iterate once
    // i.e. can replace the loop with a constant variable definition
    bool dummy_loop_;
    
  };
  
};

#endif

