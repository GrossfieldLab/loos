debug_opts='-g -Wall'
release_opts='-O3'


env = Environment()
release=ARGUMENTS.get('release', 0)
if int(release):
    env.Append(CCFLAGS=release_opts)
else:
    env.Append(CCFLAGS=debug_opts)

Export('env')

SConscript('SConscript')

