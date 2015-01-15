module.exports = function (grunt) {

	grunt.initConfig({
		pkg: grunt.file.readJSON('package.json'),
		build_dir: 'dist',
		build_name: 'TagSelect',

		/* SASS COMPILATION
		**********************************/
		sass: {
			options: { style: 'compressed' },
			main: {
				files: {'<%= build_dir %>/<%= build_name %>.min.css': 'src/<%= build_name %>.scss'}
			}
		},

		/* JS MINIFICATION
		**********************************/
		uglify: {
			main: {
				options: { banner: '/*! <%= pkg.name %> <%= grunt.template.today("yyyy-mm-dd") %> */\n' },
				files: {
					'<%= build_dir %>/<%= build_name %>.min.js': 'src/*.js'
				}
			}
		},

		/* SEMVER HELPER
		**********************************/
		bump: {
			options: {
				pushTo: 'origin',
				files: ['package.json', 'bower.json'],
				commitFiles: ['package.json', 'bower.json']
			}
		},

		/* WATCH
		**********************************/
		watch: {
			sass: {
				files: 'src/*.scss',
				tasks: ['sass']
			}
		}

	});

	grunt.loadNpmTasks('grunt-contrib-uglify');
	grunt.loadNpmTasks('grunt-contrib-watch');
	grunt.loadNpmTasks('grunt-contrib-sass');
	grunt.loadNpmTasks('grunt-bump');

	grunt.registerTask('default', ['sass']);
	grunt.registerTask('build', ['default', 'uglify']);
};
