# TagSelect.js

HTML5-based input for typing lists of tags, tokens, or other discrete values.

> **Compatibility:** Only modern browsers are supported. Some of the nicer features of this plugin (copy/paste support, native selection and cursor movement) aren't realistic on older browsers.

***

## Quickstart

First, include the files from `dist/` in your project. They'll need to be included *after* jQuery. Example:

```html
<link rel="stylesheet" href="TagSelect.min.css">
<script src="TagSelect.min.js"></script>
```

Next, you can use the plugin in your own JavaScript. If you have an element with the ID 'my-input':

```html
<div id="my-input"></div>
```

```javascript
var tags = new TagSelect({
    $el: $('#my-input'),
    options: ['Red', 'Green', 'Blue', 'Yellow', 'Orange', 'Magenta']
});
```

### Methods

```javascript
/* .get()
*************************/
var selected = tags.get();
// => ['Magenta, 'Yellow']

/* .set()
*************************/
tags.set(['Orange']);

/* .on()
*************************/
tags.on('change', function () {
    console.log('Tags updated!');
    console.log(tags.get());
})
```

## Contributing

### Development

To get started, you'll need install the development dependencies:

```
npm install
```

Development tasks, like compiling [Sass](http://sass-lang.com/) and minifying JavaScript, are handled by [Grunt](http://gruntjs.com/).

```bash
# Watch for changes to *.scss files and recompile
grunt watch

# Rebuild minified CSS and JS files
grunt build
```
