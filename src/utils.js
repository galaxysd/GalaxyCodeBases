var utils;

if (process.env.BUILD === 'lite') {
  utils = {
    escape: _.escape,
    invert: _.invert,
    isEqual: _.isEqual,
    merge: _.merge,
    uniq: _.uniq
  };
} else {
  utils = {
    escape: require('escape-html'),
    invert: require('array-to-map'),
    isEqual: require('array-equal'),
    merge: require('object-assign'),
    uniq: require('array-uniq')
  };
}

utils.escapeStringRegexp = require('escape-string-regexp');

/**
 * Thanks to Tim Down.
 * http://stackoverflow.com/a/4812022/1314762
 *
 * Linebreaks and some CSS are not handled, so if
 *  that becomes an issue we'll likely need to use
 *  the larger (45kb) rangy-core.js plugin.
 */
utils.getCaretCharacterOffsetWithin = function (element) {
  var caretOffset = 0;
  var doc = element.ownerDocument || element.document;
  var win = doc.defaultView || doc.parentWindow;
  var sel;
  if (typeof win.getSelection !== 'undefined') {
    sel = win.getSelection();
    if (sel.rangeCount > 0) {
      var range = win.getSelection().getRangeAt(0);
      var preCaretRange = range.cloneRange();
      preCaretRange.selectNodeContents(element);
      preCaretRange.setEnd(range.endContainer, range.endOffset);
      caretOffset = preCaretRange.toString().length;
    }
  } else if ( (sel = doc.selection) && sel.type !== 'Control') {
    var textRange = sel.createRange();
    var preCaretTextRange = doc.body.createTextRange();
    preCaretTextRange.moveToElementText(element);
    preCaretTextRange.setEndPoint('EndToEnd', textRange);
    caretOffset = preCaretTextRange.text.length;
  }
  return caretOffset;
};

module.exports = utils;
