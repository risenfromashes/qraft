#include "internal.hpp"

#include <cctype>

namespace qraft {
namespace {

/// Newick parser. The internal nodes still open are kept on an explicit
/// stack, so the nesting depth of the tree costs no call stack.
class Parser {
public:
  explicit Parser(std::string_view text) : text_(text) {}

  /// Parses exactly one tree. The closing ';' is optional; only whitespace
  /// and comments may follow it.
  Tree parse() {
    Tree tree;
    // Every node except the root follows a '(' or a ',', so this bounds the
    // node count (exactly, unless labels or comments contain those
    // characters) and parsed trees hold no growth slack.
    std::size_t node_bound = 1;
    for (char ch : text_) node_bound += ch == '(' || ch == ',' ? 1 : 0;
    tree.nodes.reserve(node_bound);
    skip_ws_and_comments();
    const int root = parse_subtree(tree);
    tree.root = root;
    skip_ws_and_comments();
    if (peek() == ';') {
      ++pos_;
    }
    skip_ws_and_comments();
    if (pos_ != text_.size()) {
      throw std::invalid_argument("unexpected trailing Newick text");
    }
    return tree;
  }

private:
  std::string_view text_;
  std::size_t pos_ = 0;

  /// The current character, or '\0' at the end of the text.
  char peek() const { return pos_ < text_.size() ? text_[pos_] : '\0'; }

  /// Consumes the current character; throws at the end of the text.
  char get() {
    if (pos_ >= text_.size()) {
      throw std::invalid_argument("unexpected end of Newick text");
    }
    return text_[pos_++];
  }

  /// Skips whitespace and bracketed comments, which may nest.
  void skip_ws_and_comments() {
    bool progressed = true;
    while (progressed) {
      progressed = false;
      while (std::isspace(static_cast<unsigned char>(peek()))) {
        ++pos_;
        progressed = true;
      }
      if (peek() == '[') {
        ++pos_;
        int depth = 1;
        while (pos_ < text_.size() && depth > 0) {
          const char ch = text_[pos_++];
          if (ch == '[') {
            ++depth;
          } else if (ch == ']') {
            --depth;
          }
        }
        if (depth != 0) {
          throw std::invalid_argument("unterminated Newick comment");
        }
        progressed = true;
      }
    }
  }

  int add_node(Tree &tree, std::string label = {}) {
    const int id = static_cast<int>(tree.nodes.size());
    tree.nodes.push_back(Tree::Node{{}, std::move(label)});
    return id;
  }

  void add_edge(Tree &tree, int a, int b) {
    tree.nodes[a].neighbors.push_back(b);
    tree.nodes[b].neighbors.push_back(a);
  }

  /// Parses a subtree and returns its root. An internal node is created at
  /// its '(' and joined to each child once the child is complete; its label
  /// and branch length are read and discarded. A leaf must have a label.
  int parse_subtree(Tree &tree) {
    std::vector<int> open;
    while (true) {
      // Start a subtree: open an internal node, or read a leaf.
      skip_ws_and_comments();
      if (peek() == '(') {
        get();
        open.push_back(add_node(tree));
        continue;
      }
      const std::string label = parse_label();
      if (label.empty()) {
        throw std::invalid_argument("leaf without label in Newick tree");
      }
      skip_branch_length();
      int completed = add_node(tree, label);

      // Join the completed subtree to its parent, closing every parent
      // whose ')' follows, until a ',' starts a sibling.
      while (true) {
        if (open.empty()) {
          return completed;
        }
        const int node = open.back();
        add_edge(tree, node, completed);
        skip_ws_and_comments();
        if (peek() == ',') {
          get();
          break;
        }
        if (peek() != ')') {
          throw std::invalid_argument(
              "expected ',' or ')' in Newick subtree");
        }
        get();
        parse_label();
        skip_branch_length();
        open.pop_back();
        completed = node;
      }
    }
  }

  /// Parses a quoted or unquoted label. Inside quotes, '' stands for one
  /// quote character.
  std::string parse_label() {
    skip_ws_and_comments();
    if (peek() == '\'') {
      get();
      std::string label;
      while (true) {
        const char ch = get();
        if (ch == '\'') {
          if (peek() == '\'') {
            get();
            label += '\'';
            continue;
          }
          break;
        }
        label += ch;
      }
      skip_ws_and_comments();
      return label;
    }

    std::string label;
    while (pos_ < text_.size()) {
      const char ch = peek();
      if (ch == ':' || ch == ',' || ch == ')' || ch == '(' || ch == ';' ||
          ch == '[' || ch == ']' ||
          std::isspace(static_cast<unsigned char>(ch))) {
        break;
      }
      label += ch;
      ++pos_;
    }
    skip_ws_and_comments();
    return label;
  }

  /// Skips an optional ':' and the branch length after it.
  void skip_branch_length() {
    skip_ws_and_comments();
    if (peek() != ':') {
      return;
    }
    ++pos_;
    while (pos_ < text_.size()) {
      const char ch = peek();
      if (ch == ',' || ch == ')' || ch == ';' || ch == '[') {
        break;
      }
      ++pos_;
    }
    skip_ws_and_comments();
  }
};

} // namespace

Tree parse_newick(std::string_view text) { return Parser(text).parse(); }

} // namespace qraft
