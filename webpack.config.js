var webpack = require('webpack'),
	LodashModuleReplacementPlugin = require('lodash-webpack-plugin');

process.env.BUILD = (process.argv.indexOf('--lite') !== -1) ? 'lite' : 'full';

var plugins = [
	new webpack.EnvironmentPlugin(['BUILD']),
	new LodashModuleReplacementPlugin(),
	new webpack.optimize.OccurenceOrderPlugin()
];

if (process.argv.indexOf('--minify') !== -1) {
	plugins.push(new webpack.optimize.UglifyJsPlugin());
}

module.exports = {
	'entry': './browser.js',
	'module': {
		'loaders': [{
			'loader': 'babel',
			'test': /\.js$/,
			'exclude': /node_modules/,
			'query': {
				'plugins': ['lodash'],
				'presets': ['es2015']
			}
		}],
		'postLoaders': [{
			loader: 'transform?envify'
		}]
	},
	'plugins': plugins
};
