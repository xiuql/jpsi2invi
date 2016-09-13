
#-- start of constituents_header ------

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

tags      = $(tag),$(CMTEXTRATAGS)

TestRelease_tag = $(tag)

#cmt_local_tagfile = $(TestRelease_tag).make
cmt_local_tagfile = $(bin)$(TestRelease_tag).make

#-include $(cmt_local_tagfile)
include $(cmt_local_tagfile)

#cmt_local_setup = $(bin)setup$$$$.make
#cmt_local_setup = $(bin)$(package)setup$$$$.make
#cmt_final_setup = $(bin)TestReleasesetup.make
cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)$(package)setup.make

cmt_build_library_linksstamp = $(bin)cmt_build_library_links.stamp
#--------------------------------------------------------

#cmt_lock_setup = /tmp/lock$(cmt_lock_pid).make
#cmt_temp_tag = /tmp/tag$(cmt_lock_pid).make

#first :: $(cmt_local_tagfile)
#	@echo $(cmt_local_tagfile) ok
#ifndef QUICK
#first :: $(cmt_final_setup) ;
#else
#first :: ;
#endif

##	@bin=`$(cmtexe) show macro_value bin`

#$(cmt_local_tagfile) : $(cmt_lock_setup)
#	@echo "#CMT> Error: $@: No such file" >&2; exit 1
#$(cmt_local_tagfile) :
#	@echo "#CMT> Warning: $@: No such file" >&2; exit
#	@echo "#CMT> Info: $@: No need to rebuild file" >&2; exit

#$(cmt_final_setup) : $(cmt_local_tagfile) 
#	$(echo) "(constituents.make) Rebuilding $@"
#	@if test ! -d $(@D); then $(mkdir) -p $(@D); fi; \
#	  if test -f $(cmt_local_setup); then /bin/rm -f $(cmt_local_setup); fi; \
#	  trap '/bin/rm -f $(cmt_local_setup)' 0 1 2 15; \
#	  $(cmtexe) -tag=$(tags) show setup >>$(cmt_local_setup); \
#	  if test ! -f $@; then \
#	    mv $(cmt_local_setup) $@; \
#	  else \
#	    if /usr/bin/diff $(cmt_local_setup) $@ >/dev/null ; then \
#	      : ; \
#	    else \
#	      mv $(cmt_local_setup) $@; \
#	    fi; \
#	  fi

#	@/bin/echo $@ ok   

#config :: checkuses
#	@exit 0
#checkuses : ;

env.make ::
	printenv >env.make.tmp; $(cmtexe) check files env.make.tmp env.make

ifndef QUICK
all :: build_library_links ;
else
all :: $(cmt_build_library_linksstamp) ;
endif

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

dirs :: requirements
	@if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi
#	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
#	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

#requirements :
#	@if test ! -r requirements ; then echo "No requirements file" ; fi

build_library_links : dirs
	$(echo) "(constituents.make) Rebuilding library links"; \
	 $(build_library_links)
#	if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi; \
#	$(build_library_links)

$(cmt_build_library_linksstamp) : $(cmt_final_setup) $(cmt_local_tagfile) $(bin)library_links.in
	$(echo) "(constituents.make) Rebuilding library links"; \
	 $(build_library_links) -f=$(bin)library_links.in -without_cmt
	$(silent) \touch $@

ifndef PEDANTIC
.DEFAULT ::
#.DEFAULT :
	$(echo) "(constituents.make) $@: No rule for such target" >&2
endif

${CMTROOT}/src/Makefile.core : ;
ifdef use_requirements
$(use_requirements) : ;
endif

#-- end of constituents_header ------
#-- start of group ------

all_groups :: all

all :: $(all_dependencies)  $(all_pre_constituents) $(all_constituents)  $(all_post_constituents)
	$(echo) "all ok."

#	@/bin/echo " all ok."

clean :: allclean

allclean ::  $(all_constituentsclean)
	$(echo) $(all_constituentsclean)
	$(echo) "allclean ok."

#	@echo $(all_constituentsclean)
#	@/bin/echo " allclean ok."

#-- end of group ------
#-- start of group ------

all_groups :: cmt_actions

cmt_actions :: $(cmt_actions_dependencies)  $(cmt_actions_pre_constituents) $(cmt_actions_constituents)  $(cmt_actions_post_constituents)
	$(echo) "cmt_actions ok."

#	@/bin/echo " cmt_actions ok."

clean :: allclean

cmt_actionsclean ::  $(cmt_actions_constituentsclean)
	$(echo) $(cmt_actions_constituentsclean)
	$(echo) "cmt_actionsclean ok."

#	@echo $(cmt_actions_constituentsclean)
#	@/bin/echo " cmt_actionsclean ok."

#-- end of group ------
#-- start of constituent ------

cmt_config_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_config_has_target_tag

#cmt_local_tagfile_config = $(TestRelease_tag)_config.make
cmt_local_tagfile_config = $(bin)$(TestRelease_tag)_config.make
cmt_local_setup_config = $(bin)setup_config$$$$.make
cmt_final_setup_config = $(bin)setup_config.make
#cmt_final_setup_config = $(bin)TestRelease_configsetup.make
cmt_local_config_makefile = $(bin)config.make

config_extratags = -tag_add=target_config

#$(cmt_local_tagfile_config) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_config) ::
else
$(cmt_local_tagfile_config) :
endif
	$(echo) "(constituents.make) Rebuilding $@"; \
	  if test -f $(cmt_local_tagfile_config); then /bin/rm -f $(cmt_local_tagfile_config); fi ; \
	  $(cmtexe) -tag=$(tags) $(config_extratags) build tag_makefile >>$(cmt_local_tagfile_config)
	$(echo) "(constituents.make) Rebuilding $(cmt_final_setup_config)"; \
	  test ! -f $(cmt_local_setup_config) || \rm -f $(cmt_local_setup_config); \
	  trap '\rm -f $(cmt_local_setup_config)' 0 1 2 15; \
	  $(cmtexe) -tag=$(tags) $(config_extratags) show setup >$(cmt_local_setup_config) && \
	  if [ -f $(cmt_final_setup_config) ] && \
	    \cmp -s $(cmt_final_setup_config) $(cmt_local_setup_config); then \
	    \rm $(cmt_local_setup_config); else \
	    \mv -f $(cmt_local_setup_config) $(cmt_final_setup_config); fi

else

#cmt_local_tagfile_config = $(TestRelease_tag).make
cmt_local_tagfile_config = $(bin)$(TestRelease_tag).make
cmt_final_setup_config = $(bin)setup.make
#cmt_final_setup_config = $(bin)TestReleasesetup.make
cmt_local_config_makefile = $(bin)config.make

endif

not_config_dependencies = { n=0; for p in $?; do m=0; for d in $(config_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
configdirs :
	@if test ! -d $(bin)config; then $(mkdir) -p $(bin)config; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)config
else
configdirs : ;
endif

#ifndef QUICK
#ifdef STRUCTURED_OUTPUT
# configdirs ::
#	@if test ! -d $(bin)config; then $(mkdir) -p $(bin)config; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)config
#
#$(cmt_local_config_makefile) :: $(config_dependencies) $(cmt_local_tagfile_config) build_library_links dirs configdirs
#else
#$(cmt_local_config_makefile) :: $(config_dependencies) $(cmt_local_tagfile_config) build_library_links dirs
#endif
#else
#$(cmt_local_config_makefile) :: $(cmt_local_tagfile_config)
#endif

ifdef cmt_config_has_target_tag

ifndef QUICK
$(cmt_local_config_makefile) : $(config_dependencies) build_library_links
	$(echo) "(constituents.make) Building config.make"; \
	  $(cmtexe) -tag=$(tags) $(config_extratags) build constituent_config -out=$(cmt_local_config_makefile) config
else
$(cmt_local_config_makefile) : $(config_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_config) ] || \
	  [ ! -f $(cmt_final_setup_config) ] || \
	  $(not_config_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building config.make"; \
	  $(cmtexe) -tag=$(tags) $(config_extratags) build constituent_config -out=$(cmt_local_config_makefile) config; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_config_makefile) : $(config_dependencies) build_library_links
	$(echo) "(constituents.make) Building config.make"; \
	  $(cmtexe) -f=$(bin)config.in -tag=$(tags) $(config_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_config_makefile) config
else
$(cmt_local_config_makefile) : $(config_dependencies) $(cmt_build_library_linksstamp) $(bin)config.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_config) ] || \
	  [ ! -f $(cmt_final_setup_config) ] || \
	  $(not_config_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building config.make"; \
	  $(cmtexe) -f=$(bin)config.in -tag=$(tags) $(config_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_config_makefile) config; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(config_extratags) build constituent_makefile -out=$(cmt_local_config_makefile) config

config :: $(config_dependencies) $(cmt_local_config_makefile) dirs configdirs
	$(echo) "(constituents.make) Starting config"
	@if test -f $(cmt_local_config_makefile); then \
	  $(MAKE) -f $(cmt_local_config_makefile) config; \
	  fi
#	@$(MAKE) -f $(cmt_local_config_makefile) config
	$(echo) "(constituents.make) config done"

clean :: configclean

configclean :: $(configclean_dependencies) ##$(cmt_local_config_makefile)
	$(echo) "(constituents.make) Starting configclean"
	@-if test -f $(cmt_local_config_makefile); then \
	  $(MAKE) -f $(cmt_local_config_makefile) configclean; \
	fi
	$(echo) "(constituents.make) configclean done"
#	@-$(MAKE) -f $(cmt_local_config_makefile) configclean

##	  /bin/rm -f $(cmt_local_config_makefile) $(bin)config_dependencies.make

install :: configinstall

configinstall :: $(config_dependencies) $(cmt_local_config_makefile)
	$(echo) "(constituents.make) Starting install config"
	@-$(MAKE) -f $(cmt_local_config_makefile) install
	$(echo) "(constituents.make) install config done"

uninstall : configuninstall

$(foreach d,$(config_dependencies),$(eval $(d)uninstall_dependencies += configuninstall))

configuninstall : $(configuninstall_dependencies) ##$(cmt_local_config_makefile)
	$(echo) "(constituents.make) Starting uninstall config"
	@if test -f $(cmt_local_config_makefile); then \
	  $(MAKE) -f $(cmt_local_config_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_config_makefile) uninstall
	$(echo) "(constituents.make) uninstall config done"

remove_library_links :: configuninstall

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ config"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ config done"
endif

#-- end of constituent ------
#-- start of constituent ------

cmt_install_includes_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_install_includes_has_target_tag

#cmt_local_tagfile_install_includes = $(TestRelease_tag)_install_includes.make
cmt_local_tagfile_install_includes = $(bin)$(TestRelease_tag)_install_includes.make
cmt_local_setup_install_includes = $(bin)setup_install_includes$$$$.make
cmt_final_setup_install_includes = $(bin)setup_install_includes.make
#cmt_final_setup_install_includes = $(bin)TestRelease_install_includessetup.make
cmt_local_install_includes_makefile = $(bin)install_includes.make

install_includes_extratags = -tag_add=target_install_includes

#$(cmt_local_tagfile_install_includes) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_install_includes) ::
else
$(cmt_local_tagfile_install_includes) :
endif
	$(echo) "(constituents.make) Rebuilding $@"; \
	  if test -f $(cmt_local_tagfile_install_includes); then /bin/rm -f $(cmt_local_tagfile_install_includes); fi ; \
	  $(cmtexe) -tag=$(tags) $(install_includes_extratags) build tag_makefile >>$(cmt_local_tagfile_install_includes)
	$(echo) "(constituents.make) Rebuilding $(cmt_final_setup_install_includes)"; \
	  test ! -f $(cmt_local_setup_install_includes) || \rm -f $(cmt_local_setup_install_includes); \
	  trap '\rm -f $(cmt_local_setup_install_includes)' 0 1 2 15; \
	  $(cmtexe) -tag=$(tags) $(install_includes_extratags) show setup >$(cmt_local_setup_install_includes) && \
	  if [ -f $(cmt_final_setup_install_includes) ] && \
	    \cmp -s $(cmt_final_setup_install_includes) $(cmt_local_setup_install_includes); then \
	    \rm $(cmt_local_setup_install_includes); else \
	    \mv -f $(cmt_local_setup_install_includes) $(cmt_final_setup_install_includes); fi

else

#cmt_local_tagfile_install_includes = $(TestRelease_tag).make
cmt_local_tagfile_install_includes = $(bin)$(TestRelease_tag).make
cmt_final_setup_install_includes = $(bin)setup.make
#cmt_final_setup_install_includes = $(bin)TestReleasesetup.make
cmt_local_install_includes_makefile = $(bin)install_includes.make

endif

not_install_includes_dependencies = { n=0; for p in $?; do m=0; for d in $(install_includes_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
install_includesdirs :
	@if test ! -d $(bin)install_includes; then $(mkdir) -p $(bin)install_includes; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)install_includes
else
install_includesdirs : ;
endif

#ifndef QUICK
#ifdef STRUCTURED_OUTPUT
# install_includesdirs ::
#	@if test ! -d $(bin)install_includes; then $(mkdir) -p $(bin)install_includes; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)install_includes
#
#$(cmt_local_install_includes_makefile) :: $(install_includes_dependencies) $(cmt_local_tagfile_install_includes) build_library_links dirs install_includesdirs
#else
#$(cmt_local_install_includes_makefile) :: $(install_includes_dependencies) $(cmt_local_tagfile_install_includes) build_library_links dirs
#endif
#else
#$(cmt_local_install_includes_makefile) :: $(cmt_local_tagfile_install_includes)
#endif

ifdef cmt_install_includes_has_target_tag

ifndef QUICK
$(cmt_local_install_includes_makefile) : $(install_includes_dependencies) build_library_links
	$(echo) "(constituents.make) Building install_includes.make"; \
	  $(cmtexe) -tag=$(tags) $(install_includes_extratags) build constituent_config -out=$(cmt_local_install_includes_makefile) install_includes
else
$(cmt_local_install_includes_makefile) : $(install_includes_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_install_includes) ] || \
	  [ ! -f $(cmt_final_setup_install_includes) ] || \
	  $(not_install_includes_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building install_includes.make"; \
	  $(cmtexe) -tag=$(tags) $(install_includes_extratags) build constituent_config -out=$(cmt_local_install_includes_makefile) install_includes; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_install_includes_makefile) : $(install_includes_dependencies) build_library_links
	$(echo) "(constituents.make) Building install_includes.make"; \
	  $(cmtexe) -f=$(bin)install_includes.in -tag=$(tags) $(install_includes_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_install_includes_makefile) install_includes
else
$(cmt_local_install_includes_makefile) : $(install_includes_dependencies) $(cmt_build_library_linksstamp) $(bin)install_includes.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_install_includes) ] || \
	  [ ! -f $(cmt_final_setup_install_includes) ] || \
	  $(not_install_includes_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building install_includes.make"; \
	  $(cmtexe) -f=$(bin)install_includes.in -tag=$(tags) $(install_includes_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_install_includes_makefile) install_includes; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(install_includes_extratags) build constituent_makefile -out=$(cmt_local_install_includes_makefile) install_includes

install_includes :: $(install_includes_dependencies) $(cmt_local_install_includes_makefile) dirs install_includesdirs
	$(echo) "(constituents.make) Starting install_includes"
	@if test -f $(cmt_local_install_includes_makefile); then \
	  $(MAKE) -f $(cmt_local_install_includes_makefile) install_includes; \
	  fi
#	@$(MAKE) -f $(cmt_local_install_includes_makefile) install_includes
	$(echo) "(constituents.make) install_includes done"

clean :: install_includesclean

install_includesclean :: $(install_includesclean_dependencies) ##$(cmt_local_install_includes_makefile)
	$(echo) "(constituents.make) Starting install_includesclean"
	@-if test -f $(cmt_local_install_includes_makefile); then \
	  $(MAKE) -f $(cmt_local_install_includes_makefile) install_includesclean; \
	fi
	$(echo) "(constituents.make) install_includesclean done"
#	@-$(MAKE) -f $(cmt_local_install_includes_makefile) install_includesclean

##	  /bin/rm -f $(cmt_local_install_includes_makefile) $(bin)install_includes_dependencies.make

install :: install_includesinstall

install_includesinstall :: $(install_includes_dependencies) $(cmt_local_install_includes_makefile)
	$(echo) "(constituents.make) Starting install install_includes"
	@-$(MAKE) -f $(cmt_local_install_includes_makefile) install
	$(echo) "(constituents.make) install install_includes done"

uninstall : install_includesuninstall

$(foreach d,$(install_includes_dependencies),$(eval $(d)uninstall_dependencies += install_includesuninstall))

install_includesuninstall : $(install_includesuninstall_dependencies) ##$(cmt_local_install_includes_makefile)
	$(echo) "(constituents.make) Starting uninstall install_includes"
	@if test -f $(cmt_local_install_includes_makefile); then \
	  $(MAKE) -f $(cmt_local_install_includes_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_install_includes_makefile) uninstall
	$(echo) "(constituents.make) uninstall install_includes done"

remove_library_links :: install_includesuninstall

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ install_includes"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ install_includes done"
endif

#-- end of constituent ------
#-- start of constituent ------

cmt_make_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_make_has_target_tag

#cmt_local_tagfile_make = $(TestRelease_tag)_make.make
cmt_local_tagfile_make = $(bin)$(TestRelease_tag)_make.make
cmt_local_setup_make = $(bin)setup_make$$$$.make
cmt_final_setup_make = $(bin)setup_make.make
#cmt_final_setup_make = $(bin)TestRelease_makesetup.make
cmt_local_make_makefile = $(bin)make.make

make_extratags = -tag_add=target_make

#$(cmt_local_tagfile_make) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_make) ::
else
$(cmt_local_tagfile_make) :
endif
	$(echo) "(constituents.make) Rebuilding $@"; \
	  if test -f $(cmt_local_tagfile_make); then /bin/rm -f $(cmt_local_tagfile_make); fi ; \
	  $(cmtexe) -tag=$(tags) $(make_extratags) build tag_makefile >>$(cmt_local_tagfile_make)
	$(echo) "(constituents.make) Rebuilding $(cmt_final_setup_make)"; \
	  test ! -f $(cmt_local_setup_make) || \rm -f $(cmt_local_setup_make); \
	  trap '\rm -f $(cmt_local_setup_make)' 0 1 2 15; \
	  $(cmtexe) -tag=$(tags) $(make_extratags) show setup >$(cmt_local_setup_make) && \
	  if [ -f $(cmt_final_setup_make) ] && \
	    \cmp -s $(cmt_final_setup_make) $(cmt_local_setup_make); then \
	    \rm $(cmt_local_setup_make); else \
	    \mv -f $(cmt_local_setup_make) $(cmt_final_setup_make); fi

else

#cmt_local_tagfile_make = $(TestRelease_tag).make
cmt_local_tagfile_make = $(bin)$(TestRelease_tag).make
cmt_final_setup_make = $(bin)setup.make
#cmt_final_setup_make = $(bin)TestReleasesetup.make
cmt_local_make_makefile = $(bin)make.make

endif

not_make_dependencies = { n=0; for p in $?; do m=0; for d in $(make_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
makedirs :
	@if test ! -d $(bin)make; then $(mkdir) -p $(bin)make; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)make
else
makedirs : ;
endif

#ifndef QUICK
#ifdef STRUCTURED_OUTPUT
# makedirs ::
#	@if test ! -d $(bin)make; then $(mkdir) -p $(bin)make; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)make
#
#$(cmt_local_make_makefile) :: $(make_dependencies) $(cmt_local_tagfile_make) build_library_links dirs makedirs
#else
#$(cmt_local_make_makefile) :: $(make_dependencies) $(cmt_local_tagfile_make) build_library_links dirs
#endif
#else
#$(cmt_local_make_makefile) :: $(cmt_local_tagfile_make)
#endif

ifdef cmt_make_has_target_tag

ifndef QUICK
$(cmt_local_make_makefile) : $(make_dependencies) build_library_links
	$(echo) "(constituents.make) Building make.make"; \
	  $(cmtexe) -tag=$(tags) $(make_extratags) build constituent_config -out=$(cmt_local_make_makefile) make
else
$(cmt_local_make_makefile) : $(make_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_make) ] || \
	  [ ! -f $(cmt_final_setup_make) ] || \
	  $(not_make_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building make.make"; \
	  $(cmtexe) -tag=$(tags) $(make_extratags) build constituent_config -out=$(cmt_local_make_makefile) make; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_make_makefile) : $(make_dependencies) build_library_links
	$(echo) "(constituents.make) Building make.make"; \
	  $(cmtexe) -f=$(bin)make.in -tag=$(tags) $(make_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_make_makefile) make
else
$(cmt_local_make_makefile) : $(make_dependencies) $(cmt_build_library_linksstamp) $(bin)make.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_make) ] || \
	  [ ! -f $(cmt_final_setup_make) ] || \
	  $(not_make_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building make.make"; \
	  $(cmtexe) -f=$(bin)make.in -tag=$(tags) $(make_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_make_makefile) make; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(make_extratags) build constituent_makefile -out=$(cmt_local_make_makefile) make

make :: $(make_dependencies) $(cmt_local_make_makefile) dirs makedirs
	$(echo) "(constituents.make) Starting make"
	@if test -f $(cmt_local_make_makefile); then \
	  $(MAKE) -f $(cmt_local_make_makefile) make; \
	  fi
#	@$(MAKE) -f $(cmt_local_make_makefile) make
	$(echo) "(constituents.make) make done"

clean :: makeclean

makeclean :: $(makeclean_dependencies) ##$(cmt_local_make_makefile)
	$(echo) "(constituents.make) Starting makeclean"
	@-if test -f $(cmt_local_make_makefile); then \
	  $(MAKE) -f $(cmt_local_make_makefile) makeclean; \
	fi
	$(echo) "(constituents.make) makeclean done"
#	@-$(MAKE) -f $(cmt_local_make_makefile) makeclean

##	  /bin/rm -f $(cmt_local_make_makefile) $(bin)make_dependencies.make

install :: makeinstall

makeinstall :: $(make_dependencies) $(cmt_local_make_makefile)
	$(echo) "(constituents.make) Starting install make"
	@-$(MAKE) -f $(cmt_local_make_makefile) install
	$(echo) "(constituents.make) install make done"

uninstall : makeuninstall

$(foreach d,$(make_dependencies),$(eval $(d)uninstall_dependencies += makeuninstall))

makeuninstall : $(makeuninstall_dependencies) ##$(cmt_local_make_makefile)
	$(echo) "(constituents.make) Starting uninstall make"
	@if test -f $(cmt_local_make_makefile); then \
	  $(MAKE) -f $(cmt_local_make_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_make_makefile) uninstall
	$(echo) "(constituents.make) uninstall make done"

remove_library_links :: makeuninstall

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ make"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ make done"
endif

#-- end of constituent ------
#-- start of constituent ------

cmt_CompilePython_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_CompilePython_has_target_tag

#cmt_local_tagfile_CompilePython = $(TestRelease_tag)_CompilePython.make
cmt_local_tagfile_CompilePython = $(bin)$(TestRelease_tag)_CompilePython.make
cmt_local_setup_CompilePython = $(bin)setup_CompilePython$$$$.make
cmt_final_setup_CompilePython = $(bin)setup_CompilePython.make
#cmt_final_setup_CompilePython = $(bin)TestRelease_CompilePythonsetup.make
cmt_local_CompilePython_makefile = $(bin)CompilePython.make

CompilePython_extratags = -tag_add=target_CompilePython

#$(cmt_local_tagfile_CompilePython) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_CompilePython) ::
else
$(cmt_local_tagfile_CompilePython) :
endif
	$(echo) "(constituents.make) Rebuilding $@"; \
	  if test -f $(cmt_local_tagfile_CompilePython); then /bin/rm -f $(cmt_local_tagfile_CompilePython); fi ; \
	  $(cmtexe) -tag=$(tags) $(CompilePython_extratags) build tag_makefile >>$(cmt_local_tagfile_CompilePython)
	$(echo) "(constituents.make) Rebuilding $(cmt_final_setup_CompilePython)"; \
	  test ! -f $(cmt_local_setup_CompilePython) || \rm -f $(cmt_local_setup_CompilePython); \
	  trap '\rm -f $(cmt_local_setup_CompilePython)' 0 1 2 15; \
	  $(cmtexe) -tag=$(tags) $(CompilePython_extratags) show setup >$(cmt_local_setup_CompilePython) && \
	  if [ -f $(cmt_final_setup_CompilePython) ] && \
	    \cmp -s $(cmt_final_setup_CompilePython) $(cmt_local_setup_CompilePython); then \
	    \rm $(cmt_local_setup_CompilePython); else \
	    \mv -f $(cmt_local_setup_CompilePython) $(cmt_final_setup_CompilePython); fi

else

#cmt_local_tagfile_CompilePython = $(TestRelease_tag).make
cmt_local_tagfile_CompilePython = $(bin)$(TestRelease_tag).make
cmt_final_setup_CompilePython = $(bin)setup.make
#cmt_final_setup_CompilePython = $(bin)TestReleasesetup.make
cmt_local_CompilePython_makefile = $(bin)CompilePython.make

endif

not_CompilePython_dependencies = { n=0; for p in $?; do m=0; for d in $(CompilePython_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
CompilePythondirs :
	@if test ! -d $(bin)CompilePython; then $(mkdir) -p $(bin)CompilePython; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)CompilePython
else
CompilePythondirs : ;
endif

#ifndef QUICK
#ifdef STRUCTURED_OUTPUT
# CompilePythondirs ::
#	@if test ! -d $(bin)CompilePython; then $(mkdir) -p $(bin)CompilePython; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)CompilePython
#
#$(cmt_local_CompilePython_makefile) :: $(CompilePython_dependencies) $(cmt_local_tagfile_CompilePython) build_library_links dirs CompilePythondirs
#else
#$(cmt_local_CompilePython_makefile) :: $(CompilePython_dependencies) $(cmt_local_tagfile_CompilePython) build_library_links dirs
#endif
#else
#$(cmt_local_CompilePython_makefile) :: $(cmt_local_tagfile_CompilePython)
#endif

ifdef cmt_CompilePython_has_target_tag

ifndef QUICK
$(cmt_local_CompilePython_makefile) : $(CompilePython_dependencies) build_library_links
	$(echo) "(constituents.make) Building CompilePython.make"; \
	  $(cmtexe) -tag=$(tags) $(CompilePython_extratags) build constituent_config -out=$(cmt_local_CompilePython_makefile) CompilePython
else
$(cmt_local_CompilePython_makefile) : $(CompilePython_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_CompilePython) ] || \
	  [ ! -f $(cmt_final_setup_CompilePython) ] || \
	  $(not_CompilePython_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building CompilePython.make"; \
	  $(cmtexe) -tag=$(tags) $(CompilePython_extratags) build constituent_config -out=$(cmt_local_CompilePython_makefile) CompilePython; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_CompilePython_makefile) : $(CompilePython_dependencies) build_library_links
	$(echo) "(constituents.make) Building CompilePython.make"; \
	  $(cmtexe) -f=$(bin)CompilePython.in -tag=$(tags) $(CompilePython_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_CompilePython_makefile) CompilePython
else
$(cmt_local_CompilePython_makefile) : $(CompilePython_dependencies) $(cmt_build_library_linksstamp) $(bin)CompilePython.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_CompilePython) ] || \
	  [ ! -f $(cmt_final_setup_CompilePython) ] || \
	  $(not_CompilePython_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building CompilePython.make"; \
	  $(cmtexe) -f=$(bin)CompilePython.in -tag=$(tags) $(CompilePython_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_CompilePython_makefile) CompilePython; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(CompilePython_extratags) build constituent_makefile -out=$(cmt_local_CompilePython_makefile) CompilePython

CompilePython :: $(CompilePython_dependencies) $(cmt_local_CompilePython_makefile) dirs CompilePythondirs
	$(echo) "(constituents.make) Starting CompilePython"
	@if test -f $(cmt_local_CompilePython_makefile); then \
	  $(MAKE) -f $(cmt_local_CompilePython_makefile) CompilePython; \
	  fi
#	@$(MAKE) -f $(cmt_local_CompilePython_makefile) CompilePython
	$(echo) "(constituents.make) CompilePython done"

clean :: CompilePythonclean

CompilePythonclean :: $(CompilePythonclean_dependencies) ##$(cmt_local_CompilePython_makefile)
	$(echo) "(constituents.make) Starting CompilePythonclean"
	@-if test -f $(cmt_local_CompilePython_makefile); then \
	  $(MAKE) -f $(cmt_local_CompilePython_makefile) CompilePythonclean; \
	fi
	$(echo) "(constituents.make) CompilePythonclean done"
#	@-$(MAKE) -f $(cmt_local_CompilePython_makefile) CompilePythonclean

##	  /bin/rm -f $(cmt_local_CompilePython_makefile) $(bin)CompilePython_dependencies.make

install :: CompilePythoninstall

CompilePythoninstall :: $(CompilePython_dependencies) $(cmt_local_CompilePython_makefile)
	$(echo) "(constituents.make) Starting install CompilePython"
	@-$(MAKE) -f $(cmt_local_CompilePython_makefile) install
	$(echo) "(constituents.make) install CompilePython done"

uninstall : CompilePythonuninstall

$(foreach d,$(CompilePython_dependencies),$(eval $(d)uninstall_dependencies += CompilePythonuninstall))

CompilePythonuninstall : $(CompilePythonuninstall_dependencies) ##$(cmt_local_CompilePython_makefile)
	$(echo) "(constituents.make) Starting uninstall CompilePython"
	@if test -f $(cmt_local_CompilePython_makefile); then \
	  $(MAKE) -f $(cmt_local_CompilePython_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_CompilePython_makefile) uninstall
	$(echo) "(constituents.make) uninstall CompilePython done"

remove_library_links :: CompilePythonuninstall

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ CompilePython"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ CompilePython done"
endif

#-- end of constituent ------
#-- start of constituent ------

cmt_qmtest_run_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_qmtest_run_has_target_tag

#cmt_local_tagfile_qmtest_run = $(TestRelease_tag)_qmtest_run.make
cmt_local_tagfile_qmtest_run = $(bin)$(TestRelease_tag)_qmtest_run.make
cmt_local_setup_qmtest_run = $(bin)setup_qmtest_run$$$$.make
cmt_final_setup_qmtest_run = $(bin)setup_qmtest_run.make
#cmt_final_setup_qmtest_run = $(bin)TestRelease_qmtest_runsetup.make
cmt_local_qmtest_run_makefile = $(bin)qmtest_run.make

qmtest_run_extratags = -tag_add=target_qmtest_run

#$(cmt_local_tagfile_qmtest_run) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_qmtest_run) ::
else
$(cmt_local_tagfile_qmtest_run) :
endif
	$(echo) "(constituents.make) Rebuilding $@"; \
	  if test -f $(cmt_local_tagfile_qmtest_run); then /bin/rm -f $(cmt_local_tagfile_qmtest_run); fi ; \
	  $(cmtexe) -tag=$(tags) $(qmtest_run_extratags) build tag_makefile >>$(cmt_local_tagfile_qmtest_run)
	$(echo) "(constituents.make) Rebuilding $(cmt_final_setup_qmtest_run)"; \
	  test ! -f $(cmt_local_setup_qmtest_run) || \rm -f $(cmt_local_setup_qmtest_run); \
	  trap '\rm -f $(cmt_local_setup_qmtest_run)' 0 1 2 15; \
	  $(cmtexe) -tag=$(tags) $(qmtest_run_extratags) show setup >$(cmt_local_setup_qmtest_run) && \
	  if [ -f $(cmt_final_setup_qmtest_run) ] && \
	    \cmp -s $(cmt_final_setup_qmtest_run) $(cmt_local_setup_qmtest_run); then \
	    \rm $(cmt_local_setup_qmtest_run); else \
	    \mv -f $(cmt_local_setup_qmtest_run) $(cmt_final_setup_qmtest_run); fi

else

#cmt_local_tagfile_qmtest_run = $(TestRelease_tag).make
cmt_local_tagfile_qmtest_run = $(bin)$(TestRelease_tag).make
cmt_final_setup_qmtest_run = $(bin)setup.make
#cmt_final_setup_qmtest_run = $(bin)TestReleasesetup.make
cmt_local_qmtest_run_makefile = $(bin)qmtest_run.make

endif

not_qmtest_run_dependencies = { n=0; for p in $?; do m=0; for d in $(qmtest_run_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
qmtest_rundirs :
	@if test ! -d $(bin)qmtest_run; then $(mkdir) -p $(bin)qmtest_run; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)qmtest_run
else
qmtest_rundirs : ;
endif

#ifndef QUICK
#ifdef STRUCTURED_OUTPUT
# qmtest_rundirs ::
#	@if test ! -d $(bin)qmtest_run; then $(mkdir) -p $(bin)qmtest_run; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)qmtest_run
#
#$(cmt_local_qmtest_run_makefile) :: $(qmtest_run_dependencies) $(cmt_local_tagfile_qmtest_run) build_library_links dirs qmtest_rundirs
#else
#$(cmt_local_qmtest_run_makefile) :: $(qmtest_run_dependencies) $(cmt_local_tagfile_qmtest_run) build_library_links dirs
#endif
#else
#$(cmt_local_qmtest_run_makefile) :: $(cmt_local_tagfile_qmtest_run)
#endif

ifdef cmt_qmtest_run_has_target_tag

ifndef QUICK
$(cmt_local_qmtest_run_makefile) : $(qmtest_run_dependencies) build_library_links
	$(echo) "(constituents.make) Building qmtest_run.make"; \
	  $(cmtexe) -tag=$(tags) $(qmtest_run_extratags) build constituent_config -out=$(cmt_local_qmtest_run_makefile) qmtest_run
else
$(cmt_local_qmtest_run_makefile) : $(qmtest_run_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_qmtest_run) ] || \
	  [ ! -f $(cmt_final_setup_qmtest_run) ] || \
	  $(not_qmtest_run_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building qmtest_run.make"; \
	  $(cmtexe) -tag=$(tags) $(qmtest_run_extratags) build constituent_config -out=$(cmt_local_qmtest_run_makefile) qmtest_run; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_qmtest_run_makefile) : $(qmtest_run_dependencies) build_library_links
	$(echo) "(constituents.make) Building qmtest_run.make"; \
	  $(cmtexe) -f=$(bin)qmtest_run.in -tag=$(tags) $(qmtest_run_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_qmtest_run_makefile) qmtest_run
else
$(cmt_local_qmtest_run_makefile) : $(qmtest_run_dependencies) $(cmt_build_library_linksstamp) $(bin)qmtest_run.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_qmtest_run) ] || \
	  [ ! -f $(cmt_final_setup_qmtest_run) ] || \
	  $(not_qmtest_run_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building qmtest_run.make"; \
	  $(cmtexe) -f=$(bin)qmtest_run.in -tag=$(tags) $(qmtest_run_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_qmtest_run_makefile) qmtest_run; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(qmtest_run_extratags) build constituent_makefile -out=$(cmt_local_qmtest_run_makefile) qmtest_run

qmtest_run :: $(qmtest_run_dependencies) $(cmt_local_qmtest_run_makefile) dirs qmtest_rundirs
	$(echo) "(constituents.make) Starting qmtest_run"
	@if test -f $(cmt_local_qmtest_run_makefile); then \
	  $(MAKE) -f $(cmt_local_qmtest_run_makefile) qmtest_run; \
	  fi
#	@$(MAKE) -f $(cmt_local_qmtest_run_makefile) qmtest_run
	$(echo) "(constituents.make) qmtest_run done"

clean :: qmtest_runclean

qmtest_runclean :: $(qmtest_runclean_dependencies) ##$(cmt_local_qmtest_run_makefile)
	$(echo) "(constituents.make) Starting qmtest_runclean"
	@-if test -f $(cmt_local_qmtest_run_makefile); then \
	  $(MAKE) -f $(cmt_local_qmtest_run_makefile) qmtest_runclean; \
	fi
	$(echo) "(constituents.make) qmtest_runclean done"
#	@-$(MAKE) -f $(cmt_local_qmtest_run_makefile) qmtest_runclean

##	  /bin/rm -f $(cmt_local_qmtest_run_makefile) $(bin)qmtest_run_dependencies.make

install :: qmtest_runinstall

qmtest_runinstall :: $(qmtest_run_dependencies) $(cmt_local_qmtest_run_makefile)
	$(echo) "(constituents.make) Starting install qmtest_run"
	@-$(MAKE) -f $(cmt_local_qmtest_run_makefile) install
	$(echo) "(constituents.make) install qmtest_run done"

uninstall : qmtest_rununinstall

$(foreach d,$(qmtest_run_dependencies),$(eval $(d)uninstall_dependencies += qmtest_rununinstall))

qmtest_rununinstall : $(qmtest_rununinstall_dependencies) ##$(cmt_local_qmtest_run_makefile)
	$(echo) "(constituents.make) Starting uninstall qmtest_run"
	@if test -f $(cmt_local_qmtest_run_makefile); then \
	  $(MAKE) -f $(cmt_local_qmtest_run_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_qmtest_run_makefile) uninstall
	$(echo) "(constituents.make) uninstall qmtest_run done"

remove_library_links :: qmtest_rununinstall

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ qmtest_run"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ qmtest_run done"
endif

#-- end of constituent ------
#-- start of constituent ------

cmt_qmtest_summarize_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_qmtest_summarize_has_target_tag

#cmt_local_tagfile_qmtest_summarize = $(TestRelease_tag)_qmtest_summarize.make
cmt_local_tagfile_qmtest_summarize = $(bin)$(TestRelease_tag)_qmtest_summarize.make
cmt_local_setup_qmtest_summarize = $(bin)setup_qmtest_summarize$$$$.make
cmt_final_setup_qmtest_summarize = $(bin)setup_qmtest_summarize.make
#cmt_final_setup_qmtest_summarize = $(bin)TestRelease_qmtest_summarizesetup.make
cmt_local_qmtest_summarize_makefile = $(bin)qmtest_summarize.make

qmtest_summarize_extratags = -tag_add=target_qmtest_summarize

#$(cmt_local_tagfile_qmtest_summarize) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_qmtest_summarize) ::
else
$(cmt_local_tagfile_qmtest_summarize) :
endif
	$(echo) "(constituents.make) Rebuilding $@"; \
	  if test -f $(cmt_local_tagfile_qmtest_summarize); then /bin/rm -f $(cmt_local_tagfile_qmtest_summarize); fi ; \
	  $(cmtexe) -tag=$(tags) $(qmtest_summarize_extratags) build tag_makefile >>$(cmt_local_tagfile_qmtest_summarize)
	$(echo) "(constituents.make) Rebuilding $(cmt_final_setup_qmtest_summarize)"; \
	  test ! -f $(cmt_local_setup_qmtest_summarize) || \rm -f $(cmt_local_setup_qmtest_summarize); \
	  trap '\rm -f $(cmt_local_setup_qmtest_summarize)' 0 1 2 15; \
	  $(cmtexe) -tag=$(tags) $(qmtest_summarize_extratags) show setup >$(cmt_local_setup_qmtest_summarize) && \
	  if [ -f $(cmt_final_setup_qmtest_summarize) ] && \
	    \cmp -s $(cmt_final_setup_qmtest_summarize) $(cmt_local_setup_qmtest_summarize); then \
	    \rm $(cmt_local_setup_qmtest_summarize); else \
	    \mv -f $(cmt_local_setup_qmtest_summarize) $(cmt_final_setup_qmtest_summarize); fi

else

#cmt_local_tagfile_qmtest_summarize = $(TestRelease_tag).make
cmt_local_tagfile_qmtest_summarize = $(bin)$(TestRelease_tag).make
cmt_final_setup_qmtest_summarize = $(bin)setup.make
#cmt_final_setup_qmtest_summarize = $(bin)TestReleasesetup.make
cmt_local_qmtest_summarize_makefile = $(bin)qmtest_summarize.make

endif

not_qmtest_summarize_dependencies = { n=0; for p in $?; do m=0; for d in $(qmtest_summarize_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
qmtest_summarizedirs :
	@if test ! -d $(bin)qmtest_summarize; then $(mkdir) -p $(bin)qmtest_summarize; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)qmtest_summarize
else
qmtest_summarizedirs : ;
endif

#ifndef QUICK
#ifdef STRUCTURED_OUTPUT
# qmtest_summarizedirs ::
#	@if test ! -d $(bin)qmtest_summarize; then $(mkdir) -p $(bin)qmtest_summarize; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)qmtest_summarize
#
#$(cmt_local_qmtest_summarize_makefile) :: $(qmtest_summarize_dependencies) $(cmt_local_tagfile_qmtest_summarize) build_library_links dirs qmtest_summarizedirs
#else
#$(cmt_local_qmtest_summarize_makefile) :: $(qmtest_summarize_dependencies) $(cmt_local_tagfile_qmtest_summarize) build_library_links dirs
#endif
#else
#$(cmt_local_qmtest_summarize_makefile) :: $(cmt_local_tagfile_qmtest_summarize)
#endif

ifdef cmt_qmtest_summarize_has_target_tag

ifndef QUICK
$(cmt_local_qmtest_summarize_makefile) : $(qmtest_summarize_dependencies) build_library_links
	$(echo) "(constituents.make) Building qmtest_summarize.make"; \
	  $(cmtexe) -tag=$(tags) $(qmtest_summarize_extratags) build constituent_config -out=$(cmt_local_qmtest_summarize_makefile) qmtest_summarize
else
$(cmt_local_qmtest_summarize_makefile) : $(qmtest_summarize_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_qmtest_summarize) ] || \
	  [ ! -f $(cmt_final_setup_qmtest_summarize) ] || \
	  $(not_qmtest_summarize_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building qmtest_summarize.make"; \
	  $(cmtexe) -tag=$(tags) $(qmtest_summarize_extratags) build constituent_config -out=$(cmt_local_qmtest_summarize_makefile) qmtest_summarize; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_qmtest_summarize_makefile) : $(qmtest_summarize_dependencies) build_library_links
	$(echo) "(constituents.make) Building qmtest_summarize.make"; \
	  $(cmtexe) -f=$(bin)qmtest_summarize.in -tag=$(tags) $(qmtest_summarize_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_qmtest_summarize_makefile) qmtest_summarize
else
$(cmt_local_qmtest_summarize_makefile) : $(qmtest_summarize_dependencies) $(cmt_build_library_linksstamp) $(bin)qmtest_summarize.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_qmtest_summarize) ] || \
	  [ ! -f $(cmt_final_setup_qmtest_summarize) ] || \
	  $(not_qmtest_summarize_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building qmtest_summarize.make"; \
	  $(cmtexe) -f=$(bin)qmtest_summarize.in -tag=$(tags) $(qmtest_summarize_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_qmtest_summarize_makefile) qmtest_summarize; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(qmtest_summarize_extratags) build constituent_makefile -out=$(cmt_local_qmtest_summarize_makefile) qmtest_summarize

qmtest_summarize :: $(qmtest_summarize_dependencies) $(cmt_local_qmtest_summarize_makefile) dirs qmtest_summarizedirs
	$(echo) "(constituents.make) Starting qmtest_summarize"
	@if test -f $(cmt_local_qmtest_summarize_makefile); then \
	  $(MAKE) -f $(cmt_local_qmtest_summarize_makefile) qmtest_summarize; \
	  fi
#	@$(MAKE) -f $(cmt_local_qmtest_summarize_makefile) qmtest_summarize
	$(echo) "(constituents.make) qmtest_summarize done"

clean :: qmtest_summarizeclean

qmtest_summarizeclean :: $(qmtest_summarizeclean_dependencies) ##$(cmt_local_qmtest_summarize_makefile)
	$(echo) "(constituents.make) Starting qmtest_summarizeclean"
	@-if test -f $(cmt_local_qmtest_summarize_makefile); then \
	  $(MAKE) -f $(cmt_local_qmtest_summarize_makefile) qmtest_summarizeclean; \
	fi
	$(echo) "(constituents.make) qmtest_summarizeclean done"
#	@-$(MAKE) -f $(cmt_local_qmtest_summarize_makefile) qmtest_summarizeclean

##	  /bin/rm -f $(cmt_local_qmtest_summarize_makefile) $(bin)qmtest_summarize_dependencies.make

install :: qmtest_summarizeinstall

qmtest_summarizeinstall :: $(qmtest_summarize_dependencies) $(cmt_local_qmtest_summarize_makefile)
	$(echo) "(constituents.make) Starting install qmtest_summarize"
	@-$(MAKE) -f $(cmt_local_qmtest_summarize_makefile) install
	$(echo) "(constituents.make) install qmtest_summarize done"

uninstall : qmtest_summarizeuninstall

$(foreach d,$(qmtest_summarize_dependencies),$(eval $(d)uninstall_dependencies += qmtest_summarizeuninstall))

qmtest_summarizeuninstall : $(qmtest_summarizeuninstall_dependencies) ##$(cmt_local_qmtest_summarize_makefile)
	$(echo) "(constituents.make) Starting uninstall qmtest_summarize"
	@if test -f $(cmt_local_qmtest_summarize_makefile); then \
	  $(MAKE) -f $(cmt_local_qmtest_summarize_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_qmtest_summarize_makefile) uninstall
	$(echo) "(constituents.make) uninstall qmtest_summarize done"

remove_library_links :: qmtest_summarizeuninstall

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ qmtest_summarize"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ qmtest_summarize done"
endif

#-- end of constituent ------
#-- start of constituent ------

cmt_TestPackage_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_TestPackage_has_target_tag

#cmt_local_tagfile_TestPackage = $(TestRelease_tag)_TestPackage.make
cmt_local_tagfile_TestPackage = $(bin)$(TestRelease_tag)_TestPackage.make
cmt_local_setup_TestPackage = $(bin)setup_TestPackage$$$$.make
cmt_final_setup_TestPackage = $(bin)setup_TestPackage.make
#cmt_final_setup_TestPackage = $(bin)TestRelease_TestPackagesetup.make
cmt_local_TestPackage_makefile = $(bin)TestPackage.make

TestPackage_extratags = -tag_add=target_TestPackage

#$(cmt_local_tagfile_TestPackage) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_TestPackage) ::
else
$(cmt_local_tagfile_TestPackage) :
endif
	$(echo) "(constituents.make) Rebuilding $@"; \
	  if test -f $(cmt_local_tagfile_TestPackage); then /bin/rm -f $(cmt_local_tagfile_TestPackage); fi ; \
	  $(cmtexe) -tag=$(tags) $(TestPackage_extratags) build tag_makefile >>$(cmt_local_tagfile_TestPackage)
	$(echo) "(constituents.make) Rebuilding $(cmt_final_setup_TestPackage)"; \
	  test ! -f $(cmt_local_setup_TestPackage) || \rm -f $(cmt_local_setup_TestPackage); \
	  trap '\rm -f $(cmt_local_setup_TestPackage)' 0 1 2 15; \
	  $(cmtexe) -tag=$(tags) $(TestPackage_extratags) show setup >$(cmt_local_setup_TestPackage) && \
	  if [ -f $(cmt_final_setup_TestPackage) ] && \
	    \cmp -s $(cmt_final_setup_TestPackage) $(cmt_local_setup_TestPackage); then \
	    \rm $(cmt_local_setup_TestPackage); else \
	    \mv -f $(cmt_local_setup_TestPackage) $(cmt_final_setup_TestPackage); fi

else

#cmt_local_tagfile_TestPackage = $(TestRelease_tag).make
cmt_local_tagfile_TestPackage = $(bin)$(TestRelease_tag).make
cmt_final_setup_TestPackage = $(bin)setup.make
#cmt_final_setup_TestPackage = $(bin)TestReleasesetup.make
cmt_local_TestPackage_makefile = $(bin)TestPackage.make

endif

not_TestPackage_dependencies = { n=0; for p in $?; do m=0; for d in $(TestPackage_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
TestPackagedirs :
	@if test ! -d $(bin)TestPackage; then $(mkdir) -p $(bin)TestPackage; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)TestPackage
else
TestPackagedirs : ;
endif

#ifndef QUICK
#ifdef STRUCTURED_OUTPUT
# TestPackagedirs ::
#	@if test ! -d $(bin)TestPackage; then $(mkdir) -p $(bin)TestPackage; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)TestPackage
#
#$(cmt_local_TestPackage_makefile) :: $(TestPackage_dependencies) $(cmt_local_tagfile_TestPackage) build_library_links dirs TestPackagedirs
#else
#$(cmt_local_TestPackage_makefile) :: $(TestPackage_dependencies) $(cmt_local_tagfile_TestPackage) build_library_links dirs
#endif
#else
#$(cmt_local_TestPackage_makefile) :: $(cmt_local_tagfile_TestPackage)
#endif

ifdef cmt_TestPackage_has_target_tag

ifndef QUICK
$(cmt_local_TestPackage_makefile) : $(TestPackage_dependencies) build_library_links
	$(echo) "(constituents.make) Building TestPackage.make"; \
	  $(cmtexe) -tag=$(tags) $(TestPackage_extratags) build constituent_config -out=$(cmt_local_TestPackage_makefile) TestPackage
else
$(cmt_local_TestPackage_makefile) : $(TestPackage_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_TestPackage) ] || \
	  [ ! -f $(cmt_final_setup_TestPackage) ] || \
	  $(not_TestPackage_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building TestPackage.make"; \
	  $(cmtexe) -tag=$(tags) $(TestPackage_extratags) build constituent_config -out=$(cmt_local_TestPackage_makefile) TestPackage; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_TestPackage_makefile) : $(TestPackage_dependencies) build_library_links
	$(echo) "(constituents.make) Building TestPackage.make"; \
	  $(cmtexe) -f=$(bin)TestPackage.in -tag=$(tags) $(TestPackage_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_TestPackage_makefile) TestPackage
else
$(cmt_local_TestPackage_makefile) : $(TestPackage_dependencies) $(cmt_build_library_linksstamp) $(bin)TestPackage.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_TestPackage) ] || \
	  [ ! -f $(cmt_final_setup_TestPackage) ] || \
	  $(not_TestPackage_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building TestPackage.make"; \
	  $(cmtexe) -f=$(bin)TestPackage.in -tag=$(tags) $(TestPackage_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_TestPackage_makefile) TestPackage; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(TestPackage_extratags) build constituent_makefile -out=$(cmt_local_TestPackage_makefile) TestPackage

TestPackage :: $(TestPackage_dependencies) $(cmt_local_TestPackage_makefile) dirs TestPackagedirs
	$(echo) "(constituents.make) Starting TestPackage"
	@if test -f $(cmt_local_TestPackage_makefile); then \
	  $(MAKE) -f $(cmt_local_TestPackage_makefile) TestPackage; \
	  fi
#	@$(MAKE) -f $(cmt_local_TestPackage_makefile) TestPackage
	$(echo) "(constituents.make) TestPackage done"

clean :: TestPackageclean

TestPackageclean :: $(TestPackageclean_dependencies) ##$(cmt_local_TestPackage_makefile)
	$(echo) "(constituents.make) Starting TestPackageclean"
	@-if test -f $(cmt_local_TestPackage_makefile); then \
	  $(MAKE) -f $(cmt_local_TestPackage_makefile) TestPackageclean; \
	fi
	$(echo) "(constituents.make) TestPackageclean done"
#	@-$(MAKE) -f $(cmt_local_TestPackage_makefile) TestPackageclean

##	  /bin/rm -f $(cmt_local_TestPackage_makefile) $(bin)TestPackage_dependencies.make

install :: TestPackageinstall

TestPackageinstall :: $(TestPackage_dependencies) $(cmt_local_TestPackage_makefile)
	$(echo) "(constituents.make) Starting install TestPackage"
	@-$(MAKE) -f $(cmt_local_TestPackage_makefile) install
	$(echo) "(constituents.make) install TestPackage done"

uninstall : TestPackageuninstall

$(foreach d,$(TestPackage_dependencies),$(eval $(d)uninstall_dependencies += TestPackageuninstall))

TestPackageuninstall : $(TestPackageuninstall_dependencies) ##$(cmt_local_TestPackage_makefile)
	$(echo) "(constituents.make) Starting uninstall TestPackage"
	@if test -f $(cmt_local_TestPackage_makefile); then \
	  $(MAKE) -f $(cmt_local_TestPackage_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_TestPackage_makefile) uninstall
	$(echo) "(constituents.make) uninstall TestPackage done"

remove_library_links :: TestPackageuninstall

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ TestPackage"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ TestPackage done"
endif

#-- end of constituent ------
#-- start of constituent ------

cmt_TestProject_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_TestProject_has_target_tag

#cmt_local_tagfile_TestProject = $(TestRelease_tag)_TestProject.make
cmt_local_tagfile_TestProject = $(bin)$(TestRelease_tag)_TestProject.make
cmt_local_setup_TestProject = $(bin)setup_TestProject$$$$.make
cmt_final_setup_TestProject = $(bin)setup_TestProject.make
#cmt_final_setup_TestProject = $(bin)TestRelease_TestProjectsetup.make
cmt_local_TestProject_makefile = $(bin)TestProject.make

TestProject_extratags = -tag_add=target_TestProject

#$(cmt_local_tagfile_TestProject) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_TestProject) ::
else
$(cmt_local_tagfile_TestProject) :
endif
	$(echo) "(constituents.make) Rebuilding $@"; \
	  if test -f $(cmt_local_tagfile_TestProject); then /bin/rm -f $(cmt_local_tagfile_TestProject); fi ; \
	  $(cmtexe) -tag=$(tags) $(TestProject_extratags) build tag_makefile >>$(cmt_local_tagfile_TestProject)
	$(echo) "(constituents.make) Rebuilding $(cmt_final_setup_TestProject)"; \
	  test ! -f $(cmt_local_setup_TestProject) || \rm -f $(cmt_local_setup_TestProject); \
	  trap '\rm -f $(cmt_local_setup_TestProject)' 0 1 2 15; \
	  $(cmtexe) -tag=$(tags) $(TestProject_extratags) show setup >$(cmt_local_setup_TestProject) && \
	  if [ -f $(cmt_final_setup_TestProject) ] && \
	    \cmp -s $(cmt_final_setup_TestProject) $(cmt_local_setup_TestProject); then \
	    \rm $(cmt_local_setup_TestProject); else \
	    \mv -f $(cmt_local_setup_TestProject) $(cmt_final_setup_TestProject); fi

else

#cmt_local_tagfile_TestProject = $(TestRelease_tag).make
cmt_local_tagfile_TestProject = $(bin)$(TestRelease_tag).make
cmt_final_setup_TestProject = $(bin)setup.make
#cmt_final_setup_TestProject = $(bin)TestReleasesetup.make
cmt_local_TestProject_makefile = $(bin)TestProject.make

endif

not_TestProject_dependencies = { n=0; for p in $?; do m=0; for d in $(TestProject_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
TestProjectdirs :
	@if test ! -d $(bin)TestProject; then $(mkdir) -p $(bin)TestProject; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)TestProject
else
TestProjectdirs : ;
endif

#ifndef QUICK
#ifdef STRUCTURED_OUTPUT
# TestProjectdirs ::
#	@if test ! -d $(bin)TestProject; then $(mkdir) -p $(bin)TestProject; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)TestProject
#
#$(cmt_local_TestProject_makefile) :: $(TestProject_dependencies) $(cmt_local_tagfile_TestProject) build_library_links dirs TestProjectdirs
#else
#$(cmt_local_TestProject_makefile) :: $(TestProject_dependencies) $(cmt_local_tagfile_TestProject) build_library_links dirs
#endif
#else
#$(cmt_local_TestProject_makefile) :: $(cmt_local_tagfile_TestProject)
#endif

ifdef cmt_TestProject_has_target_tag

ifndef QUICK
$(cmt_local_TestProject_makefile) : $(TestProject_dependencies) build_library_links
	$(echo) "(constituents.make) Building TestProject.make"; \
	  $(cmtexe) -tag=$(tags) $(TestProject_extratags) build constituent_config -out=$(cmt_local_TestProject_makefile) TestProject
else
$(cmt_local_TestProject_makefile) : $(TestProject_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_TestProject) ] || \
	  [ ! -f $(cmt_final_setup_TestProject) ] || \
	  $(not_TestProject_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building TestProject.make"; \
	  $(cmtexe) -tag=$(tags) $(TestProject_extratags) build constituent_config -out=$(cmt_local_TestProject_makefile) TestProject; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_TestProject_makefile) : $(TestProject_dependencies) build_library_links
	$(echo) "(constituents.make) Building TestProject.make"; \
	  $(cmtexe) -f=$(bin)TestProject.in -tag=$(tags) $(TestProject_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_TestProject_makefile) TestProject
else
$(cmt_local_TestProject_makefile) : $(TestProject_dependencies) $(cmt_build_library_linksstamp) $(bin)TestProject.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_TestProject) ] || \
	  [ ! -f $(cmt_final_setup_TestProject) ] || \
	  $(not_TestProject_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building TestProject.make"; \
	  $(cmtexe) -f=$(bin)TestProject.in -tag=$(tags) $(TestProject_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_TestProject_makefile) TestProject; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(TestProject_extratags) build constituent_makefile -out=$(cmt_local_TestProject_makefile) TestProject

TestProject :: $(TestProject_dependencies) $(cmt_local_TestProject_makefile) dirs TestProjectdirs
	$(echo) "(constituents.make) Starting TestProject"
	@if test -f $(cmt_local_TestProject_makefile); then \
	  $(MAKE) -f $(cmt_local_TestProject_makefile) TestProject; \
	  fi
#	@$(MAKE) -f $(cmt_local_TestProject_makefile) TestProject
	$(echo) "(constituents.make) TestProject done"

clean :: TestProjectclean

TestProjectclean :: $(TestProjectclean_dependencies) ##$(cmt_local_TestProject_makefile)
	$(echo) "(constituents.make) Starting TestProjectclean"
	@-if test -f $(cmt_local_TestProject_makefile); then \
	  $(MAKE) -f $(cmt_local_TestProject_makefile) TestProjectclean; \
	fi
	$(echo) "(constituents.make) TestProjectclean done"
#	@-$(MAKE) -f $(cmt_local_TestProject_makefile) TestProjectclean

##	  /bin/rm -f $(cmt_local_TestProject_makefile) $(bin)TestProject_dependencies.make

install :: TestProjectinstall

TestProjectinstall :: $(TestProject_dependencies) $(cmt_local_TestProject_makefile)
	$(echo) "(constituents.make) Starting install TestProject"
	@-$(MAKE) -f $(cmt_local_TestProject_makefile) install
	$(echo) "(constituents.make) install TestProject done"

uninstall : TestProjectuninstall

$(foreach d,$(TestProject_dependencies),$(eval $(d)uninstall_dependencies += TestProjectuninstall))

TestProjectuninstall : $(TestProjectuninstall_dependencies) ##$(cmt_local_TestProject_makefile)
	$(echo) "(constituents.make) Starting uninstall TestProject"
	@if test -f $(cmt_local_TestProject_makefile); then \
	  $(MAKE) -f $(cmt_local_TestProject_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_TestProject_makefile) uninstall
	$(echo) "(constituents.make) uninstall TestProject done"

remove_library_links :: TestProjectuninstall

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ TestProject"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ TestProject done"
endif

#-- end of constituent ------
#-- start of constituent ------

cmt_new_rootsys_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_new_rootsys_has_target_tag

#cmt_local_tagfile_new_rootsys = $(TestRelease_tag)_new_rootsys.make
cmt_local_tagfile_new_rootsys = $(bin)$(TestRelease_tag)_new_rootsys.make
cmt_local_setup_new_rootsys = $(bin)setup_new_rootsys$$$$.make
cmt_final_setup_new_rootsys = $(bin)setup_new_rootsys.make
#cmt_final_setup_new_rootsys = $(bin)TestRelease_new_rootsyssetup.make
cmt_local_new_rootsys_makefile = $(bin)new_rootsys.make

new_rootsys_extratags = -tag_add=target_new_rootsys

#$(cmt_local_tagfile_new_rootsys) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_new_rootsys) ::
else
$(cmt_local_tagfile_new_rootsys) :
endif
	$(echo) "(constituents.make) Rebuilding $@"; \
	  if test -f $(cmt_local_tagfile_new_rootsys); then /bin/rm -f $(cmt_local_tagfile_new_rootsys); fi ; \
	  $(cmtexe) -tag=$(tags) $(new_rootsys_extratags) build tag_makefile >>$(cmt_local_tagfile_new_rootsys)
	$(echo) "(constituents.make) Rebuilding $(cmt_final_setup_new_rootsys)"; \
	  test ! -f $(cmt_local_setup_new_rootsys) || \rm -f $(cmt_local_setup_new_rootsys); \
	  trap '\rm -f $(cmt_local_setup_new_rootsys)' 0 1 2 15; \
	  $(cmtexe) -tag=$(tags) $(new_rootsys_extratags) show setup >$(cmt_local_setup_new_rootsys) && \
	  if [ -f $(cmt_final_setup_new_rootsys) ] && \
	    \cmp -s $(cmt_final_setup_new_rootsys) $(cmt_local_setup_new_rootsys); then \
	    \rm $(cmt_local_setup_new_rootsys); else \
	    \mv -f $(cmt_local_setup_new_rootsys) $(cmt_final_setup_new_rootsys); fi

else

#cmt_local_tagfile_new_rootsys = $(TestRelease_tag).make
cmt_local_tagfile_new_rootsys = $(bin)$(TestRelease_tag).make
cmt_final_setup_new_rootsys = $(bin)setup.make
#cmt_final_setup_new_rootsys = $(bin)TestReleasesetup.make
cmt_local_new_rootsys_makefile = $(bin)new_rootsys.make

endif

not_new_rootsys_dependencies = { n=0; for p in $?; do m=0; for d in $(new_rootsys_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
new_rootsysdirs :
	@if test ! -d $(bin)new_rootsys; then $(mkdir) -p $(bin)new_rootsys; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)new_rootsys
else
new_rootsysdirs : ;
endif

#ifndef QUICK
#ifdef STRUCTURED_OUTPUT
# new_rootsysdirs ::
#	@if test ! -d $(bin)new_rootsys; then $(mkdir) -p $(bin)new_rootsys; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)new_rootsys
#
#$(cmt_local_new_rootsys_makefile) :: $(new_rootsys_dependencies) $(cmt_local_tagfile_new_rootsys) build_library_links dirs new_rootsysdirs
#else
#$(cmt_local_new_rootsys_makefile) :: $(new_rootsys_dependencies) $(cmt_local_tagfile_new_rootsys) build_library_links dirs
#endif
#else
#$(cmt_local_new_rootsys_makefile) :: $(cmt_local_tagfile_new_rootsys)
#endif

ifdef cmt_new_rootsys_has_target_tag

ifndef QUICK
$(cmt_local_new_rootsys_makefile) : $(new_rootsys_dependencies) build_library_links
	$(echo) "(constituents.make) Building new_rootsys.make"; \
	  $(cmtexe) -tag=$(tags) $(new_rootsys_extratags) build constituent_config -out=$(cmt_local_new_rootsys_makefile) new_rootsys
else
$(cmt_local_new_rootsys_makefile) : $(new_rootsys_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_new_rootsys) ] || \
	  [ ! -f $(cmt_final_setup_new_rootsys) ] || \
	  $(not_new_rootsys_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building new_rootsys.make"; \
	  $(cmtexe) -tag=$(tags) $(new_rootsys_extratags) build constituent_config -out=$(cmt_local_new_rootsys_makefile) new_rootsys; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_new_rootsys_makefile) : $(new_rootsys_dependencies) build_library_links
	$(echo) "(constituents.make) Building new_rootsys.make"; \
	  $(cmtexe) -f=$(bin)new_rootsys.in -tag=$(tags) $(new_rootsys_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_new_rootsys_makefile) new_rootsys
else
$(cmt_local_new_rootsys_makefile) : $(new_rootsys_dependencies) $(cmt_build_library_linksstamp) $(bin)new_rootsys.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_new_rootsys) ] || \
	  [ ! -f $(cmt_final_setup_new_rootsys) ] || \
	  $(not_new_rootsys_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building new_rootsys.make"; \
	  $(cmtexe) -f=$(bin)new_rootsys.in -tag=$(tags) $(new_rootsys_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_new_rootsys_makefile) new_rootsys; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(new_rootsys_extratags) build constituent_makefile -out=$(cmt_local_new_rootsys_makefile) new_rootsys

new_rootsys :: $(new_rootsys_dependencies) $(cmt_local_new_rootsys_makefile) dirs new_rootsysdirs
	$(echo) "(constituents.make) Starting new_rootsys"
	@if test -f $(cmt_local_new_rootsys_makefile); then \
	  $(MAKE) -f $(cmt_local_new_rootsys_makefile) new_rootsys; \
	  fi
#	@$(MAKE) -f $(cmt_local_new_rootsys_makefile) new_rootsys
	$(echo) "(constituents.make) new_rootsys done"

clean :: new_rootsysclean

new_rootsysclean :: $(new_rootsysclean_dependencies) ##$(cmt_local_new_rootsys_makefile)
	$(echo) "(constituents.make) Starting new_rootsysclean"
	@-if test -f $(cmt_local_new_rootsys_makefile); then \
	  $(MAKE) -f $(cmt_local_new_rootsys_makefile) new_rootsysclean; \
	fi
	$(echo) "(constituents.make) new_rootsysclean done"
#	@-$(MAKE) -f $(cmt_local_new_rootsys_makefile) new_rootsysclean

##	  /bin/rm -f $(cmt_local_new_rootsys_makefile) $(bin)new_rootsys_dependencies.make

install :: new_rootsysinstall

new_rootsysinstall :: $(new_rootsys_dependencies) $(cmt_local_new_rootsys_makefile)
	$(echo) "(constituents.make) Starting install new_rootsys"
	@-$(MAKE) -f $(cmt_local_new_rootsys_makefile) install
	$(echo) "(constituents.make) install new_rootsys done"

uninstall : new_rootsysuninstall

$(foreach d,$(new_rootsys_dependencies),$(eval $(d)uninstall_dependencies += new_rootsysuninstall))

new_rootsysuninstall : $(new_rootsysuninstall_dependencies) ##$(cmt_local_new_rootsys_makefile)
	$(echo) "(constituents.make) Starting uninstall new_rootsys"
	@if test -f $(cmt_local_new_rootsys_makefile); then \
	  $(MAKE) -f $(cmt_local_new_rootsys_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_new_rootsys_makefile) uninstall
	$(echo) "(constituents.make) uninstall new_rootsys done"

remove_library_links :: new_rootsysuninstall

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ new_rootsys"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ new_rootsys done"
endif

#-- end of constituent ------
#-- start of constituent ------

cmt_DatabaseSvc_check_install_scripts_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_DatabaseSvc_check_install_scripts_has_target_tag

#cmt_local_tagfile_DatabaseSvc_check_install_scripts = $(TestRelease_tag)_DatabaseSvc_check_install_scripts.make
cmt_local_tagfile_DatabaseSvc_check_install_scripts = $(bin)$(TestRelease_tag)_DatabaseSvc_check_install_scripts.make
cmt_local_setup_DatabaseSvc_check_install_scripts = $(bin)setup_DatabaseSvc_check_install_scripts$$$$.make
cmt_final_setup_DatabaseSvc_check_install_scripts = $(bin)setup_DatabaseSvc_check_install_scripts.make
#cmt_final_setup_DatabaseSvc_check_install_scripts = $(bin)TestRelease_DatabaseSvc_check_install_scriptssetup.make
cmt_local_DatabaseSvc_check_install_scripts_makefile = $(bin)DatabaseSvc_check_install_scripts.make

DatabaseSvc_check_install_scripts_extratags = -tag_add=target_DatabaseSvc_check_install_scripts

#$(cmt_local_tagfile_DatabaseSvc_check_install_scripts) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_DatabaseSvc_check_install_scripts) ::
else
$(cmt_local_tagfile_DatabaseSvc_check_install_scripts) :
endif
	$(echo) "(constituents.make) Rebuilding $@"; \
	  if test -f $(cmt_local_tagfile_DatabaseSvc_check_install_scripts); then /bin/rm -f $(cmt_local_tagfile_DatabaseSvc_check_install_scripts); fi ; \
	  $(cmtexe) -tag=$(tags) $(DatabaseSvc_check_install_scripts_extratags) build tag_makefile >>$(cmt_local_tagfile_DatabaseSvc_check_install_scripts)
	$(echo) "(constituents.make) Rebuilding $(cmt_final_setup_DatabaseSvc_check_install_scripts)"; \
	  test ! -f $(cmt_local_setup_DatabaseSvc_check_install_scripts) || \rm -f $(cmt_local_setup_DatabaseSvc_check_install_scripts); \
	  trap '\rm -f $(cmt_local_setup_DatabaseSvc_check_install_scripts)' 0 1 2 15; \
	  $(cmtexe) -tag=$(tags) $(DatabaseSvc_check_install_scripts_extratags) show setup >$(cmt_local_setup_DatabaseSvc_check_install_scripts) && \
	  if [ -f $(cmt_final_setup_DatabaseSvc_check_install_scripts) ] && \
	    \cmp -s $(cmt_final_setup_DatabaseSvc_check_install_scripts) $(cmt_local_setup_DatabaseSvc_check_install_scripts); then \
	    \rm $(cmt_local_setup_DatabaseSvc_check_install_scripts); else \
	    \mv -f $(cmt_local_setup_DatabaseSvc_check_install_scripts) $(cmt_final_setup_DatabaseSvc_check_install_scripts); fi

else

#cmt_local_tagfile_DatabaseSvc_check_install_scripts = $(TestRelease_tag).make
cmt_local_tagfile_DatabaseSvc_check_install_scripts = $(bin)$(TestRelease_tag).make
cmt_final_setup_DatabaseSvc_check_install_scripts = $(bin)setup.make
#cmt_final_setup_DatabaseSvc_check_install_scripts = $(bin)TestReleasesetup.make
cmt_local_DatabaseSvc_check_install_scripts_makefile = $(bin)DatabaseSvc_check_install_scripts.make

endif

not_DatabaseSvc_check_install_scripts_dependencies = { n=0; for p in $?; do m=0; for d in $(DatabaseSvc_check_install_scripts_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
DatabaseSvc_check_install_scriptsdirs :
	@if test ! -d $(bin)DatabaseSvc_check_install_scripts; then $(mkdir) -p $(bin)DatabaseSvc_check_install_scripts; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)DatabaseSvc_check_install_scripts
else
DatabaseSvc_check_install_scriptsdirs : ;
endif

#ifndef QUICK
#ifdef STRUCTURED_OUTPUT
# DatabaseSvc_check_install_scriptsdirs ::
#	@if test ! -d $(bin)DatabaseSvc_check_install_scripts; then $(mkdir) -p $(bin)DatabaseSvc_check_install_scripts; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)DatabaseSvc_check_install_scripts
#
#$(cmt_local_DatabaseSvc_check_install_scripts_makefile) :: $(DatabaseSvc_check_install_scripts_dependencies) $(cmt_local_tagfile_DatabaseSvc_check_install_scripts) build_library_links dirs DatabaseSvc_check_install_scriptsdirs
#else
#$(cmt_local_DatabaseSvc_check_install_scripts_makefile) :: $(DatabaseSvc_check_install_scripts_dependencies) $(cmt_local_tagfile_DatabaseSvc_check_install_scripts) build_library_links dirs
#endif
#else
#$(cmt_local_DatabaseSvc_check_install_scripts_makefile) :: $(cmt_local_tagfile_DatabaseSvc_check_install_scripts)
#endif

ifdef cmt_DatabaseSvc_check_install_scripts_has_target_tag

ifndef QUICK
$(cmt_local_DatabaseSvc_check_install_scripts_makefile) : $(DatabaseSvc_check_install_scripts_dependencies) build_library_links
	$(echo) "(constituents.make) Building DatabaseSvc_check_install_scripts.make"; \
	  $(cmtexe) -tag=$(tags) $(DatabaseSvc_check_install_scripts_extratags) build constituent_config -out=$(cmt_local_DatabaseSvc_check_install_scripts_makefile) DatabaseSvc_check_install_scripts
else
$(cmt_local_DatabaseSvc_check_install_scripts_makefile) : $(DatabaseSvc_check_install_scripts_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_DatabaseSvc_check_install_scripts) ] || \
	  [ ! -f $(cmt_final_setup_DatabaseSvc_check_install_scripts) ] || \
	  $(not_DatabaseSvc_check_install_scripts_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building DatabaseSvc_check_install_scripts.make"; \
	  $(cmtexe) -tag=$(tags) $(DatabaseSvc_check_install_scripts_extratags) build constituent_config -out=$(cmt_local_DatabaseSvc_check_install_scripts_makefile) DatabaseSvc_check_install_scripts; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_DatabaseSvc_check_install_scripts_makefile) : $(DatabaseSvc_check_install_scripts_dependencies) build_library_links
	$(echo) "(constituents.make) Building DatabaseSvc_check_install_scripts.make"; \
	  $(cmtexe) -f=$(bin)DatabaseSvc_check_install_scripts.in -tag=$(tags) $(DatabaseSvc_check_install_scripts_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_DatabaseSvc_check_install_scripts_makefile) DatabaseSvc_check_install_scripts
else
$(cmt_local_DatabaseSvc_check_install_scripts_makefile) : $(DatabaseSvc_check_install_scripts_dependencies) $(cmt_build_library_linksstamp) $(bin)DatabaseSvc_check_install_scripts.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_DatabaseSvc_check_install_scripts) ] || \
	  [ ! -f $(cmt_final_setup_DatabaseSvc_check_install_scripts) ] || \
	  $(not_DatabaseSvc_check_install_scripts_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building DatabaseSvc_check_install_scripts.make"; \
	  $(cmtexe) -f=$(bin)DatabaseSvc_check_install_scripts.in -tag=$(tags) $(DatabaseSvc_check_install_scripts_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_DatabaseSvc_check_install_scripts_makefile) DatabaseSvc_check_install_scripts; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(DatabaseSvc_check_install_scripts_extratags) build constituent_makefile -out=$(cmt_local_DatabaseSvc_check_install_scripts_makefile) DatabaseSvc_check_install_scripts

DatabaseSvc_check_install_scripts :: $(DatabaseSvc_check_install_scripts_dependencies) $(cmt_local_DatabaseSvc_check_install_scripts_makefile) dirs DatabaseSvc_check_install_scriptsdirs
	$(echo) "(constituents.make) Starting DatabaseSvc_check_install_scripts"
	@if test -f $(cmt_local_DatabaseSvc_check_install_scripts_makefile); then \
	  $(MAKE) -f $(cmt_local_DatabaseSvc_check_install_scripts_makefile) DatabaseSvc_check_install_scripts; \
	  fi
#	@$(MAKE) -f $(cmt_local_DatabaseSvc_check_install_scripts_makefile) DatabaseSvc_check_install_scripts
	$(echo) "(constituents.make) DatabaseSvc_check_install_scripts done"

clean :: DatabaseSvc_check_install_scriptsclean

DatabaseSvc_check_install_scriptsclean :: $(DatabaseSvc_check_install_scriptsclean_dependencies) ##$(cmt_local_DatabaseSvc_check_install_scripts_makefile)
	$(echo) "(constituents.make) Starting DatabaseSvc_check_install_scriptsclean"
	@-if test -f $(cmt_local_DatabaseSvc_check_install_scripts_makefile); then \
	  $(MAKE) -f $(cmt_local_DatabaseSvc_check_install_scripts_makefile) DatabaseSvc_check_install_scriptsclean; \
	fi
	$(echo) "(constituents.make) DatabaseSvc_check_install_scriptsclean done"
#	@-$(MAKE) -f $(cmt_local_DatabaseSvc_check_install_scripts_makefile) DatabaseSvc_check_install_scriptsclean

##	  /bin/rm -f $(cmt_local_DatabaseSvc_check_install_scripts_makefile) $(bin)DatabaseSvc_check_install_scripts_dependencies.make

install :: DatabaseSvc_check_install_scriptsinstall

DatabaseSvc_check_install_scriptsinstall :: $(DatabaseSvc_check_install_scripts_dependencies) $(cmt_local_DatabaseSvc_check_install_scripts_makefile)
	$(echo) "(constituents.make) Starting install DatabaseSvc_check_install_scripts"
	@-$(MAKE) -f $(cmt_local_DatabaseSvc_check_install_scripts_makefile) install
	$(echo) "(constituents.make) install DatabaseSvc_check_install_scripts done"

uninstall : DatabaseSvc_check_install_scriptsuninstall

$(foreach d,$(DatabaseSvc_check_install_scripts_dependencies),$(eval $(d)uninstall_dependencies += DatabaseSvc_check_install_scriptsuninstall))

DatabaseSvc_check_install_scriptsuninstall : $(DatabaseSvc_check_install_scriptsuninstall_dependencies) ##$(cmt_local_DatabaseSvc_check_install_scripts_makefile)
	$(echo) "(constituents.make) Starting uninstall DatabaseSvc_check_install_scripts"
	@if test -f $(cmt_local_DatabaseSvc_check_install_scripts_makefile); then \
	  $(MAKE) -f $(cmt_local_DatabaseSvc_check_install_scripts_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_DatabaseSvc_check_install_scripts_makefile) uninstall
	$(echo) "(constituents.make) uninstall DatabaseSvc_check_install_scripts done"

remove_library_links :: DatabaseSvc_check_install_scriptsuninstall

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ DatabaseSvc_check_install_scripts"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ DatabaseSvc_check_install_scripts done"
endif

#-- end of constituent ------
#-- start of constituent ------

cmt_RhopiAlg_check_install_runtime_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_RhopiAlg_check_install_runtime_has_target_tag

#cmt_local_tagfile_RhopiAlg_check_install_runtime = $(TestRelease_tag)_RhopiAlg_check_install_runtime.make
cmt_local_tagfile_RhopiAlg_check_install_runtime = $(bin)$(TestRelease_tag)_RhopiAlg_check_install_runtime.make
cmt_local_setup_RhopiAlg_check_install_runtime = $(bin)setup_RhopiAlg_check_install_runtime$$$$.make
cmt_final_setup_RhopiAlg_check_install_runtime = $(bin)setup_RhopiAlg_check_install_runtime.make
#cmt_final_setup_RhopiAlg_check_install_runtime = $(bin)TestRelease_RhopiAlg_check_install_runtimesetup.make
cmt_local_RhopiAlg_check_install_runtime_makefile = $(bin)RhopiAlg_check_install_runtime.make

RhopiAlg_check_install_runtime_extratags = -tag_add=target_RhopiAlg_check_install_runtime

#$(cmt_local_tagfile_RhopiAlg_check_install_runtime) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_RhopiAlg_check_install_runtime) ::
else
$(cmt_local_tagfile_RhopiAlg_check_install_runtime) :
endif
	$(echo) "(constituents.make) Rebuilding $@"; \
	  if test -f $(cmt_local_tagfile_RhopiAlg_check_install_runtime); then /bin/rm -f $(cmt_local_tagfile_RhopiAlg_check_install_runtime); fi ; \
	  $(cmtexe) -tag=$(tags) $(RhopiAlg_check_install_runtime_extratags) build tag_makefile >>$(cmt_local_tagfile_RhopiAlg_check_install_runtime)
	$(echo) "(constituents.make) Rebuilding $(cmt_final_setup_RhopiAlg_check_install_runtime)"; \
	  test ! -f $(cmt_local_setup_RhopiAlg_check_install_runtime) || \rm -f $(cmt_local_setup_RhopiAlg_check_install_runtime); \
	  trap '\rm -f $(cmt_local_setup_RhopiAlg_check_install_runtime)' 0 1 2 15; \
	  $(cmtexe) -tag=$(tags) $(RhopiAlg_check_install_runtime_extratags) show setup >$(cmt_local_setup_RhopiAlg_check_install_runtime) && \
	  if [ -f $(cmt_final_setup_RhopiAlg_check_install_runtime) ] && \
	    \cmp -s $(cmt_final_setup_RhopiAlg_check_install_runtime) $(cmt_local_setup_RhopiAlg_check_install_runtime); then \
	    \rm $(cmt_local_setup_RhopiAlg_check_install_runtime); else \
	    \mv -f $(cmt_local_setup_RhopiAlg_check_install_runtime) $(cmt_final_setup_RhopiAlg_check_install_runtime); fi

else

#cmt_local_tagfile_RhopiAlg_check_install_runtime = $(TestRelease_tag).make
cmt_local_tagfile_RhopiAlg_check_install_runtime = $(bin)$(TestRelease_tag).make
cmt_final_setup_RhopiAlg_check_install_runtime = $(bin)setup.make
#cmt_final_setup_RhopiAlg_check_install_runtime = $(bin)TestReleasesetup.make
cmt_local_RhopiAlg_check_install_runtime_makefile = $(bin)RhopiAlg_check_install_runtime.make

endif

not_RhopiAlg_check_install_runtime_dependencies = { n=0; for p in $?; do m=0; for d in $(RhopiAlg_check_install_runtime_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
RhopiAlg_check_install_runtimedirs :
	@if test ! -d $(bin)RhopiAlg_check_install_runtime; then $(mkdir) -p $(bin)RhopiAlg_check_install_runtime; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)RhopiAlg_check_install_runtime
else
RhopiAlg_check_install_runtimedirs : ;
endif

#ifndef QUICK
#ifdef STRUCTURED_OUTPUT
# RhopiAlg_check_install_runtimedirs ::
#	@if test ! -d $(bin)RhopiAlg_check_install_runtime; then $(mkdir) -p $(bin)RhopiAlg_check_install_runtime; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)RhopiAlg_check_install_runtime
#
#$(cmt_local_RhopiAlg_check_install_runtime_makefile) :: $(RhopiAlg_check_install_runtime_dependencies) $(cmt_local_tagfile_RhopiAlg_check_install_runtime) build_library_links dirs RhopiAlg_check_install_runtimedirs
#else
#$(cmt_local_RhopiAlg_check_install_runtime_makefile) :: $(RhopiAlg_check_install_runtime_dependencies) $(cmt_local_tagfile_RhopiAlg_check_install_runtime) build_library_links dirs
#endif
#else
#$(cmt_local_RhopiAlg_check_install_runtime_makefile) :: $(cmt_local_tagfile_RhopiAlg_check_install_runtime)
#endif

ifdef cmt_RhopiAlg_check_install_runtime_has_target_tag

ifndef QUICK
$(cmt_local_RhopiAlg_check_install_runtime_makefile) : $(RhopiAlg_check_install_runtime_dependencies) build_library_links
	$(echo) "(constituents.make) Building RhopiAlg_check_install_runtime.make"; \
	  $(cmtexe) -tag=$(tags) $(RhopiAlg_check_install_runtime_extratags) build constituent_config -out=$(cmt_local_RhopiAlg_check_install_runtime_makefile) RhopiAlg_check_install_runtime
else
$(cmt_local_RhopiAlg_check_install_runtime_makefile) : $(RhopiAlg_check_install_runtime_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_RhopiAlg_check_install_runtime) ] || \
	  [ ! -f $(cmt_final_setup_RhopiAlg_check_install_runtime) ] || \
	  $(not_RhopiAlg_check_install_runtime_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building RhopiAlg_check_install_runtime.make"; \
	  $(cmtexe) -tag=$(tags) $(RhopiAlg_check_install_runtime_extratags) build constituent_config -out=$(cmt_local_RhopiAlg_check_install_runtime_makefile) RhopiAlg_check_install_runtime; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_RhopiAlg_check_install_runtime_makefile) : $(RhopiAlg_check_install_runtime_dependencies) build_library_links
	$(echo) "(constituents.make) Building RhopiAlg_check_install_runtime.make"; \
	  $(cmtexe) -f=$(bin)RhopiAlg_check_install_runtime.in -tag=$(tags) $(RhopiAlg_check_install_runtime_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_RhopiAlg_check_install_runtime_makefile) RhopiAlg_check_install_runtime
else
$(cmt_local_RhopiAlg_check_install_runtime_makefile) : $(RhopiAlg_check_install_runtime_dependencies) $(cmt_build_library_linksstamp) $(bin)RhopiAlg_check_install_runtime.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_RhopiAlg_check_install_runtime) ] || \
	  [ ! -f $(cmt_final_setup_RhopiAlg_check_install_runtime) ] || \
	  $(not_RhopiAlg_check_install_runtime_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building RhopiAlg_check_install_runtime.make"; \
	  $(cmtexe) -f=$(bin)RhopiAlg_check_install_runtime.in -tag=$(tags) $(RhopiAlg_check_install_runtime_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_RhopiAlg_check_install_runtime_makefile) RhopiAlg_check_install_runtime; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(RhopiAlg_check_install_runtime_extratags) build constituent_makefile -out=$(cmt_local_RhopiAlg_check_install_runtime_makefile) RhopiAlg_check_install_runtime

RhopiAlg_check_install_runtime :: $(RhopiAlg_check_install_runtime_dependencies) $(cmt_local_RhopiAlg_check_install_runtime_makefile) dirs RhopiAlg_check_install_runtimedirs
	$(echo) "(constituents.make) Starting RhopiAlg_check_install_runtime"
	@if test -f $(cmt_local_RhopiAlg_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_RhopiAlg_check_install_runtime_makefile) RhopiAlg_check_install_runtime; \
	  fi
#	@$(MAKE) -f $(cmt_local_RhopiAlg_check_install_runtime_makefile) RhopiAlg_check_install_runtime
	$(echo) "(constituents.make) RhopiAlg_check_install_runtime done"

clean :: RhopiAlg_check_install_runtimeclean

RhopiAlg_check_install_runtimeclean :: $(RhopiAlg_check_install_runtimeclean_dependencies) ##$(cmt_local_RhopiAlg_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting RhopiAlg_check_install_runtimeclean"
	@-if test -f $(cmt_local_RhopiAlg_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_RhopiAlg_check_install_runtime_makefile) RhopiAlg_check_install_runtimeclean; \
	fi
	$(echo) "(constituents.make) RhopiAlg_check_install_runtimeclean done"
#	@-$(MAKE) -f $(cmt_local_RhopiAlg_check_install_runtime_makefile) RhopiAlg_check_install_runtimeclean

##	  /bin/rm -f $(cmt_local_RhopiAlg_check_install_runtime_makefile) $(bin)RhopiAlg_check_install_runtime_dependencies.make

install :: RhopiAlg_check_install_runtimeinstall

RhopiAlg_check_install_runtimeinstall :: $(RhopiAlg_check_install_runtime_dependencies) $(cmt_local_RhopiAlg_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting install RhopiAlg_check_install_runtime"
	@-$(MAKE) -f $(cmt_local_RhopiAlg_check_install_runtime_makefile) install
	$(echo) "(constituents.make) install RhopiAlg_check_install_runtime done"

uninstall : RhopiAlg_check_install_runtimeuninstall

$(foreach d,$(RhopiAlg_check_install_runtime_dependencies),$(eval $(d)uninstall_dependencies += RhopiAlg_check_install_runtimeuninstall))

RhopiAlg_check_install_runtimeuninstall : $(RhopiAlg_check_install_runtimeuninstall_dependencies) ##$(cmt_local_RhopiAlg_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting uninstall RhopiAlg_check_install_runtime"
	@if test -f $(cmt_local_RhopiAlg_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_RhopiAlg_check_install_runtime_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_RhopiAlg_check_install_runtime_makefile) uninstall
	$(echo) "(constituents.make) uninstall RhopiAlg_check_install_runtime done"

remove_library_links :: RhopiAlg_check_install_runtimeuninstall

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ RhopiAlg_check_install_runtime"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ RhopiAlg_check_install_runtime done"
endif

#-- end of constituent ------
#-- start of constituent ------

cmt_tofcalgsec_check_install_runtime_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_tofcalgsec_check_install_runtime_has_target_tag

#cmt_local_tagfile_tofcalgsec_check_install_runtime = $(TestRelease_tag)_tofcalgsec_check_install_runtime.make
cmt_local_tagfile_tofcalgsec_check_install_runtime = $(bin)$(TestRelease_tag)_tofcalgsec_check_install_runtime.make
cmt_local_setup_tofcalgsec_check_install_runtime = $(bin)setup_tofcalgsec_check_install_runtime$$$$.make
cmt_final_setup_tofcalgsec_check_install_runtime = $(bin)setup_tofcalgsec_check_install_runtime.make
#cmt_final_setup_tofcalgsec_check_install_runtime = $(bin)TestRelease_tofcalgsec_check_install_runtimesetup.make
cmt_local_tofcalgsec_check_install_runtime_makefile = $(bin)tofcalgsec_check_install_runtime.make

tofcalgsec_check_install_runtime_extratags = -tag_add=target_tofcalgsec_check_install_runtime

#$(cmt_local_tagfile_tofcalgsec_check_install_runtime) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_tofcalgsec_check_install_runtime) ::
else
$(cmt_local_tagfile_tofcalgsec_check_install_runtime) :
endif
	$(echo) "(constituents.make) Rebuilding $@"; \
	  if test -f $(cmt_local_tagfile_tofcalgsec_check_install_runtime); then /bin/rm -f $(cmt_local_tagfile_tofcalgsec_check_install_runtime); fi ; \
	  $(cmtexe) -tag=$(tags) $(tofcalgsec_check_install_runtime_extratags) build tag_makefile >>$(cmt_local_tagfile_tofcalgsec_check_install_runtime)
	$(echo) "(constituents.make) Rebuilding $(cmt_final_setup_tofcalgsec_check_install_runtime)"; \
	  test ! -f $(cmt_local_setup_tofcalgsec_check_install_runtime) || \rm -f $(cmt_local_setup_tofcalgsec_check_install_runtime); \
	  trap '\rm -f $(cmt_local_setup_tofcalgsec_check_install_runtime)' 0 1 2 15; \
	  $(cmtexe) -tag=$(tags) $(tofcalgsec_check_install_runtime_extratags) show setup >$(cmt_local_setup_tofcalgsec_check_install_runtime) && \
	  if [ -f $(cmt_final_setup_tofcalgsec_check_install_runtime) ] && \
	    \cmp -s $(cmt_final_setup_tofcalgsec_check_install_runtime) $(cmt_local_setup_tofcalgsec_check_install_runtime); then \
	    \rm $(cmt_local_setup_tofcalgsec_check_install_runtime); else \
	    \mv -f $(cmt_local_setup_tofcalgsec_check_install_runtime) $(cmt_final_setup_tofcalgsec_check_install_runtime); fi

else

#cmt_local_tagfile_tofcalgsec_check_install_runtime = $(TestRelease_tag).make
cmt_local_tagfile_tofcalgsec_check_install_runtime = $(bin)$(TestRelease_tag).make
cmt_final_setup_tofcalgsec_check_install_runtime = $(bin)setup.make
#cmt_final_setup_tofcalgsec_check_install_runtime = $(bin)TestReleasesetup.make
cmt_local_tofcalgsec_check_install_runtime_makefile = $(bin)tofcalgsec_check_install_runtime.make

endif

not_tofcalgsec_check_install_runtime_dependencies = { n=0; for p in $?; do m=0; for d in $(tofcalgsec_check_install_runtime_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
tofcalgsec_check_install_runtimedirs :
	@if test ! -d $(bin)tofcalgsec_check_install_runtime; then $(mkdir) -p $(bin)tofcalgsec_check_install_runtime; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)tofcalgsec_check_install_runtime
else
tofcalgsec_check_install_runtimedirs : ;
endif

#ifndef QUICK
#ifdef STRUCTURED_OUTPUT
# tofcalgsec_check_install_runtimedirs ::
#	@if test ! -d $(bin)tofcalgsec_check_install_runtime; then $(mkdir) -p $(bin)tofcalgsec_check_install_runtime; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)tofcalgsec_check_install_runtime
#
#$(cmt_local_tofcalgsec_check_install_runtime_makefile) :: $(tofcalgsec_check_install_runtime_dependencies) $(cmt_local_tagfile_tofcalgsec_check_install_runtime) build_library_links dirs tofcalgsec_check_install_runtimedirs
#else
#$(cmt_local_tofcalgsec_check_install_runtime_makefile) :: $(tofcalgsec_check_install_runtime_dependencies) $(cmt_local_tagfile_tofcalgsec_check_install_runtime) build_library_links dirs
#endif
#else
#$(cmt_local_tofcalgsec_check_install_runtime_makefile) :: $(cmt_local_tagfile_tofcalgsec_check_install_runtime)
#endif

ifdef cmt_tofcalgsec_check_install_runtime_has_target_tag

ifndef QUICK
$(cmt_local_tofcalgsec_check_install_runtime_makefile) : $(tofcalgsec_check_install_runtime_dependencies) build_library_links
	$(echo) "(constituents.make) Building tofcalgsec_check_install_runtime.make"; \
	  $(cmtexe) -tag=$(tags) $(tofcalgsec_check_install_runtime_extratags) build constituent_config -out=$(cmt_local_tofcalgsec_check_install_runtime_makefile) tofcalgsec_check_install_runtime
else
$(cmt_local_tofcalgsec_check_install_runtime_makefile) : $(tofcalgsec_check_install_runtime_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_tofcalgsec_check_install_runtime) ] || \
	  [ ! -f $(cmt_final_setup_tofcalgsec_check_install_runtime) ] || \
	  $(not_tofcalgsec_check_install_runtime_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building tofcalgsec_check_install_runtime.make"; \
	  $(cmtexe) -tag=$(tags) $(tofcalgsec_check_install_runtime_extratags) build constituent_config -out=$(cmt_local_tofcalgsec_check_install_runtime_makefile) tofcalgsec_check_install_runtime; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_tofcalgsec_check_install_runtime_makefile) : $(tofcalgsec_check_install_runtime_dependencies) build_library_links
	$(echo) "(constituents.make) Building tofcalgsec_check_install_runtime.make"; \
	  $(cmtexe) -f=$(bin)tofcalgsec_check_install_runtime.in -tag=$(tags) $(tofcalgsec_check_install_runtime_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_tofcalgsec_check_install_runtime_makefile) tofcalgsec_check_install_runtime
else
$(cmt_local_tofcalgsec_check_install_runtime_makefile) : $(tofcalgsec_check_install_runtime_dependencies) $(cmt_build_library_linksstamp) $(bin)tofcalgsec_check_install_runtime.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_tofcalgsec_check_install_runtime) ] || \
	  [ ! -f $(cmt_final_setup_tofcalgsec_check_install_runtime) ] || \
	  $(not_tofcalgsec_check_install_runtime_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building tofcalgsec_check_install_runtime.make"; \
	  $(cmtexe) -f=$(bin)tofcalgsec_check_install_runtime.in -tag=$(tags) $(tofcalgsec_check_install_runtime_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_tofcalgsec_check_install_runtime_makefile) tofcalgsec_check_install_runtime; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(tofcalgsec_check_install_runtime_extratags) build constituent_makefile -out=$(cmt_local_tofcalgsec_check_install_runtime_makefile) tofcalgsec_check_install_runtime

tofcalgsec_check_install_runtime :: $(tofcalgsec_check_install_runtime_dependencies) $(cmt_local_tofcalgsec_check_install_runtime_makefile) dirs tofcalgsec_check_install_runtimedirs
	$(echo) "(constituents.make) Starting tofcalgsec_check_install_runtime"
	@if test -f $(cmt_local_tofcalgsec_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_tofcalgsec_check_install_runtime_makefile) tofcalgsec_check_install_runtime; \
	  fi
#	@$(MAKE) -f $(cmt_local_tofcalgsec_check_install_runtime_makefile) tofcalgsec_check_install_runtime
	$(echo) "(constituents.make) tofcalgsec_check_install_runtime done"

clean :: tofcalgsec_check_install_runtimeclean

tofcalgsec_check_install_runtimeclean :: $(tofcalgsec_check_install_runtimeclean_dependencies) ##$(cmt_local_tofcalgsec_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting tofcalgsec_check_install_runtimeclean"
	@-if test -f $(cmt_local_tofcalgsec_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_tofcalgsec_check_install_runtime_makefile) tofcalgsec_check_install_runtimeclean; \
	fi
	$(echo) "(constituents.make) tofcalgsec_check_install_runtimeclean done"
#	@-$(MAKE) -f $(cmt_local_tofcalgsec_check_install_runtime_makefile) tofcalgsec_check_install_runtimeclean

##	  /bin/rm -f $(cmt_local_tofcalgsec_check_install_runtime_makefile) $(bin)tofcalgsec_check_install_runtime_dependencies.make

install :: tofcalgsec_check_install_runtimeinstall

tofcalgsec_check_install_runtimeinstall :: $(tofcalgsec_check_install_runtime_dependencies) $(cmt_local_tofcalgsec_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting install tofcalgsec_check_install_runtime"
	@-$(MAKE) -f $(cmt_local_tofcalgsec_check_install_runtime_makefile) install
	$(echo) "(constituents.make) install tofcalgsec_check_install_runtime done"

uninstall : tofcalgsec_check_install_runtimeuninstall

$(foreach d,$(tofcalgsec_check_install_runtime_dependencies),$(eval $(d)uninstall_dependencies += tofcalgsec_check_install_runtimeuninstall))

tofcalgsec_check_install_runtimeuninstall : $(tofcalgsec_check_install_runtimeuninstall_dependencies) ##$(cmt_local_tofcalgsec_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting uninstall tofcalgsec_check_install_runtime"
	@if test -f $(cmt_local_tofcalgsec_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_tofcalgsec_check_install_runtime_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_tofcalgsec_check_install_runtime_makefile) uninstall
	$(echo) "(constituents.make) uninstall tofcalgsec_check_install_runtime done"

remove_library_links :: tofcalgsec_check_install_runtimeuninstall

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ tofcalgsec_check_install_runtime"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ tofcalgsec_check_install_runtime done"
endif

#-- end of constituent ------
#-- start of constituent ------

cmt_OfflineEventLoopMgr_check_install_joboptions_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_OfflineEventLoopMgr_check_install_joboptions_has_target_tag

#cmt_local_tagfile_OfflineEventLoopMgr_check_install_joboptions = $(TestRelease_tag)_OfflineEventLoopMgr_check_install_joboptions.make
cmt_local_tagfile_OfflineEventLoopMgr_check_install_joboptions = $(bin)$(TestRelease_tag)_OfflineEventLoopMgr_check_install_joboptions.make
cmt_local_setup_OfflineEventLoopMgr_check_install_joboptions = $(bin)setup_OfflineEventLoopMgr_check_install_joboptions$$$$.make
cmt_final_setup_OfflineEventLoopMgr_check_install_joboptions = $(bin)setup_OfflineEventLoopMgr_check_install_joboptions.make
#cmt_final_setup_OfflineEventLoopMgr_check_install_joboptions = $(bin)TestRelease_OfflineEventLoopMgr_check_install_joboptionssetup.make
cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile = $(bin)OfflineEventLoopMgr_check_install_joboptions.make

OfflineEventLoopMgr_check_install_joboptions_extratags = -tag_add=target_OfflineEventLoopMgr_check_install_joboptions

#$(cmt_local_tagfile_OfflineEventLoopMgr_check_install_joboptions) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_OfflineEventLoopMgr_check_install_joboptions) ::
else
$(cmt_local_tagfile_OfflineEventLoopMgr_check_install_joboptions) :
endif
	$(echo) "(constituents.make) Rebuilding $@"; \
	  if test -f $(cmt_local_tagfile_OfflineEventLoopMgr_check_install_joboptions); then /bin/rm -f $(cmt_local_tagfile_OfflineEventLoopMgr_check_install_joboptions); fi ; \
	  $(cmtexe) -tag=$(tags) $(OfflineEventLoopMgr_check_install_joboptions_extratags) build tag_makefile >>$(cmt_local_tagfile_OfflineEventLoopMgr_check_install_joboptions)
	$(echo) "(constituents.make) Rebuilding $(cmt_final_setup_OfflineEventLoopMgr_check_install_joboptions)"; \
	  test ! -f $(cmt_local_setup_OfflineEventLoopMgr_check_install_joboptions) || \rm -f $(cmt_local_setup_OfflineEventLoopMgr_check_install_joboptions); \
	  trap '\rm -f $(cmt_local_setup_OfflineEventLoopMgr_check_install_joboptions)' 0 1 2 15; \
	  $(cmtexe) -tag=$(tags) $(OfflineEventLoopMgr_check_install_joboptions_extratags) show setup >$(cmt_local_setup_OfflineEventLoopMgr_check_install_joboptions) && \
	  if [ -f $(cmt_final_setup_OfflineEventLoopMgr_check_install_joboptions) ] && \
	    \cmp -s $(cmt_final_setup_OfflineEventLoopMgr_check_install_joboptions) $(cmt_local_setup_OfflineEventLoopMgr_check_install_joboptions); then \
	    \rm $(cmt_local_setup_OfflineEventLoopMgr_check_install_joboptions); else \
	    \mv -f $(cmt_local_setup_OfflineEventLoopMgr_check_install_joboptions) $(cmt_final_setup_OfflineEventLoopMgr_check_install_joboptions); fi

else

#cmt_local_tagfile_OfflineEventLoopMgr_check_install_joboptions = $(TestRelease_tag).make
cmt_local_tagfile_OfflineEventLoopMgr_check_install_joboptions = $(bin)$(TestRelease_tag).make
cmt_final_setup_OfflineEventLoopMgr_check_install_joboptions = $(bin)setup.make
#cmt_final_setup_OfflineEventLoopMgr_check_install_joboptions = $(bin)TestReleasesetup.make
cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile = $(bin)OfflineEventLoopMgr_check_install_joboptions.make

endif

not_OfflineEventLoopMgr_check_install_joboptions_dependencies = { n=0; for p in $?; do m=0; for d in $(OfflineEventLoopMgr_check_install_joboptions_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
OfflineEventLoopMgr_check_install_joboptionsdirs :
	@if test ! -d $(bin)OfflineEventLoopMgr_check_install_joboptions; then $(mkdir) -p $(bin)OfflineEventLoopMgr_check_install_joboptions; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)OfflineEventLoopMgr_check_install_joboptions
else
OfflineEventLoopMgr_check_install_joboptionsdirs : ;
endif

#ifndef QUICK
#ifdef STRUCTURED_OUTPUT
# OfflineEventLoopMgr_check_install_joboptionsdirs ::
#	@if test ! -d $(bin)OfflineEventLoopMgr_check_install_joboptions; then $(mkdir) -p $(bin)OfflineEventLoopMgr_check_install_joboptions; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)OfflineEventLoopMgr_check_install_joboptions
#
#$(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile) :: $(OfflineEventLoopMgr_check_install_joboptions_dependencies) $(cmt_local_tagfile_OfflineEventLoopMgr_check_install_joboptions) build_library_links dirs OfflineEventLoopMgr_check_install_joboptionsdirs
#else
#$(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile) :: $(OfflineEventLoopMgr_check_install_joboptions_dependencies) $(cmt_local_tagfile_OfflineEventLoopMgr_check_install_joboptions) build_library_links dirs
#endif
#else
#$(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile) :: $(cmt_local_tagfile_OfflineEventLoopMgr_check_install_joboptions)
#endif

ifdef cmt_OfflineEventLoopMgr_check_install_joboptions_has_target_tag

ifndef QUICK
$(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile) : $(OfflineEventLoopMgr_check_install_joboptions_dependencies) build_library_links
	$(echo) "(constituents.make) Building OfflineEventLoopMgr_check_install_joboptions.make"; \
	  $(cmtexe) -tag=$(tags) $(OfflineEventLoopMgr_check_install_joboptions_extratags) build constituent_config -out=$(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile) OfflineEventLoopMgr_check_install_joboptions
else
$(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile) : $(OfflineEventLoopMgr_check_install_joboptions_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_OfflineEventLoopMgr_check_install_joboptions) ] || \
	  [ ! -f $(cmt_final_setup_OfflineEventLoopMgr_check_install_joboptions) ] || \
	  $(not_OfflineEventLoopMgr_check_install_joboptions_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building OfflineEventLoopMgr_check_install_joboptions.make"; \
	  $(cmtexe) -tag=$(tags) $(OfflineEventLoopMgr_check_install_joboptions_extratags) build constituent_config -out=$(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile) OfflineEventLoopMgr_check_install_joboptions; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile) : $(OfflineEventLoopMgr_check_install_joboptions_dependencies) build_library_links
	$(echo) "(constituents.make) Building OfflineEventLoopMgr_check_install_joboptions.make"; \
	  $(cmtexe) -f=$(bin)OfflineEventLoopMgr_check_install_joboptions.in -tag=$(tags) $(OfflineEventLoopMgr_check_install_joboptions_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile) OfflineEventLoopMgr_check_install_joboptions
else
$(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile) : $(OfflineEventLoopMgr_check_install_joboptions_dependencies) $(cmt_build_library_linksstamp) $(bin)OfflineEventLoopMgr_check_install_joboptions.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_OfflineEventLoopMgr_check_install_joboptions) ] || \
	  [ ! -f $(cmt_final_setup_OfflineEventLoopMgr_check_install_joboptions) ] || \
	  $(not_OfflineEventLoopMgr_check_install_joboptions_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building OfflineEventLoopMgr_check_install_joboptions.make"; \
	  $(cmtexe) -f=$(bin)OfflineEventLoopMgr_check_install_joboptions.in -tag=$(tags) $(OfflineEventLoopMgr_check_install_joboptions_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile) OfflineEventLoopMgr_check_install_joboptions; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(OfflineEventLoopMgr_check_install_joboptions_extratags) build constituent_makefile -out=$(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile) OfflineEventLoopMgr_check_install_joboptions

OfflineEventLoopMgr_check_install_joboptions :: $(OfflineEventLoopMgr_check_install_joboptions_dependencies) $(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile) dirs OfflineEventLoopMgr_check_install_joboptionsdirs
	$(echo) "(constituents.make) Starting OfflineEventLoopMgr_check_install_joboptions"
	@if test -f $(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile); then \
	  $(MAKE) -f $(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile) OfflineEventLoopMgr_check_install_joboptions; \
	  fi
#	@$(MAKE) -f $(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile) OfflineEventLoopMgr_check_install_joboptions
	$(echo) "(constituents.make) OfflineEventLoopMgr_check_install_joboptions done"

clean :: OfflineEventLoopMgr_check_install_joboptionsclean

OfflineEventLoopMgr_check_install_joboptionsclean :: $(OfflineEventLoopMgr_check_install_joboptionsclean_dependencies) ##$(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile)
	$(echo) "(constituents.make) Starting OfflineEventLoopMgr_check_install_joboptionsclean"
	@-if test -f $(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile); then \
	  $(MAKE) -f $(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile) OfflineEventLoopMgr_check_install_joboptionsclean; \
	fi
	$(echo) "(constituents.make) OfflineEventLoopMgr_check_install_joboptionsclean done"
#	@-$(MAKE) -f $(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile) OfflineEventLoopMgr_check_install_joboptionsclean

##	  /bin/rm -f $(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile) $(bin)OfflineEventLoopMgr_check_install_joboptions_dependencies.make

install :: OfflineEventLoopMgr_check_install_joboptionsinstall

OfflineEventLoopMgr_check_install_joboptionsinstall :: $(OfflineEventLoopMgr_check_install_joboptions_dependencies) $(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile)
	$(echo) "(constituents.make) Starting install OfflineEventLoopMgr_check_install_joboptions"
	@-$(MAKE) -f $(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile) install
	$(echo) "(constituents.make) install OfflineEventLoopMgr_check_install_joboptions done"

uninstall : OfflineEventLoopMgr_check_install_joboptionsuninstall

$(foreach d,$(OfflineEventLoopMgr_check_install_joboptions_dependencies),$(eval $(d)uninstall_dependencies += OfflineEventLoopMgr_check_install_joboptionsuninstall))

OfflineEventLoopMgr_check_install_joboptionsuninstall : $(OfflineEventLoopMgr_check_install_joboptionsuninstall_dependencies) ##$(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile)
	$(echo) "(constituents.make) Starting uninstall OfflineEventLoopMgr_check_install_joboptions"
	@if test -f $(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile); then \
	  $(MAKE) -f $(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_OfflineEventLoopMgr_check_install_joboptions_makefile) uninstall
	$(echo) "(constituents.make) uninstall OfflineEventLoopMgr_check_install_joboptions done"

remove_library_links :: OfflineEventLoopMgr_check_install_joboptionsuninstall

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ OfflineEventLoopMgr_check_install_joboptions"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ OfflineEventLoopMgr_check_install_joboptions done"
endif

#-- end of constituent ------
#-- start of constituent ------

cmt_EventNavigator_check_install_runtime_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_EventNavigator_check_install_runtime_has_target_tag

#cmt_local_tagfile_EventNavigator_check_install_runtime = $(TestRelease_tag)_EventNavigator_check_install_runtime.make
cmt_local_tagfile_EventNavigator_check_install_runtime = $(bin)$(TestRelease_tag)_EventNavigator_check_install_runtime.make
cmt_local_setup_EventNavigator_check_install_runtime = $(bin)setup_EventNavigator_check_install_runtime$$$$.make
cmt_final_setup_EventNavigator_check_install_runtime = $(bin)setup_EventNavigator_check_install_runtime.make
#cmt_final_setup_EventNavigator_check_install_runtime = $(bin)TestRelease_EventNavigator_check_install_runtimesetup.make
cmt_local_EventNavigator_check_install_runtime_makefile = $(bin)EventNavigator_check_install_runtime.make

EventNavigator_check_install_runtime_extratags = -tag_add=target_EventNavigator_check_install_runtime

#$(cmt_local_tagfile_EventNavigator_check_install_runtime) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_EventNavigator_check_install_runtime) ::
else
$(cmt_local_tagfile_EventNavigator_check_install_runtime) :
endif
	$(echo) "(constituents.make) Rebuilding $@"; \
	  if test -f $(cmt_local_tagfile_EventNavigator_check_install_runtime); then /bin/rm -f $(cmt_local_tagfile_EventNavigator_check_install_runtime); fi ; \
	  $(cmtexe) -tag=$(tags) $(EventNavigator_check_install_runtime_extratags) build tag_makefile >>$(cmt_local_tagfile_EventNavigator_check_install_runtime)
	$(echo) "(constituents.make) Rebuilding $(cmt_final_setup_EventNavigator_check_install_runtime)"; \
	  test ! -f $(cmt_local_setup_EventNavigator_check_install_runtime) || \rm -f $(cmt_local_setup_EventNavigator_check_install_runtime); \
	  trap '\rm -f $(cmt_local_setup_EventNavigator_check_install_runtime)' 0 1 2 15; \
	  $(cmtexe) -tag=$(tags) $(EventNavigator_check_install_runtime_extratags) show setup >$(cmt_local_setup_EventNavigator_check_install_runtime) && \
	  if [ -f $(cmt_final_setup_EventNavigator_check_install_runtime) ] && \
	    \cmp -s $(cmt_final_setup_EventNavigator_check_install_runtime) $(cmt_local_setup_EventNavigator_check_install_runtime); then \
	    \rm $(cmt_local_setup_EventNavigator_check_install_runtime); else \
	    \mv -f $(cmt_local_setup_EventNavigator_check_install_runtime) $(cmt_final_setup_EventNavigator_check_install_runtime); fi

else

#cmt_local_tagfile_EventNavigator_check_install_runtime = $(TestRelease_tag).make
cmt_local_tagfile_EventNavigator_check_install_runtime = $(bin)$(TestRelease_tag).make
cmt_final_setup_EventNavigator_check_install_runtime = $(bin)setup.make
#cmt_final_setup_EventNavigator_check_install_runtime = $(bin)TestReleasesetup.make
cmt_local_EventNavigator_check_install_runtime_makefile = $(bin)EventNavigator_check_install_runtime.make

endif

not_EventNavigator_check_install_runtime_dependencies = { n=0; for p in $?; do m=0; for d in $(EventNavigator_check_install_runtime_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
EventNavigator_check_install_runtimedirs :
	@if test ! -d $(bin)EventNavigator_check_install_runtime; then $(mkdir) -p $(bin)EventNavigator_check_install_runtime; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)EventNavigator_check_install_runtime
else
EventNavigator_check_install_runtimedirs : ;
endif

#ifndef QUICK
#ifdef STRUCTURED_OUTPUT
# EventNavigator_check_install_runtimedirs ::
#	@if test ! -d $(bin)EventNavigator_check_install_runtime; then $(mkdir) -p $(bin)EventNavigator_check_install_runtime; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)EventNavigator_check_install_runtime
#
#$(cmt_local_EventNavigator_check_install_runtime_makefile) :: $(EventNavigator_check_install_runtime_dependencies) $(cmt_local_tagfile_EventNavigator_check_install_runtime) build_library_links dirs EventNavigator_check_install_runtimedirs
#else
#$(cmt_local_EventNavigator_check_install_runtime_makefile) :: $(EventNavigator_check_install_runtime_dependencies) $(cmt_local_tagfile_EventNavigator_check_install_runtime) build_library_links dirs
#endif
#else
#$(cmt_local_EventNavigator_check_install_runtime_makefile) :: $(cmt_local_tagfile_EventNavigator_check_install_runtime)
#endif

ifdef cmt_EventNavigator_check_install_runtime_has_target_tag

ifndef QUICK
$(cmt_local_EventNavigator_check_install_runtime_makefile) : $(EventNavigator_check_install_runtime_dependencies) build_library_links
	$(echo) "(constituents.make) Building EventNavigator_check_install_runtime.make"; \
	  $(cmtexe) -tag=$(tags) $(EventNavigator_check_install_runtime_extratags) build constituent_config -out=$(cmt_local_EventNavigator_check_install_runtime_makefile) EventNavigator_check_install_runtime
else
$(cmt_local_EventNavigator_check_install_runtime_makefile) : $(EventNavigator_check_install_runtime_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_EventNavigator_check_install_runtime) ] || \
	  [ ! -f $(cmt_final_setup_EventNavigator_check_install_runtime) ] || \
	  $(not_EventNavigator_check_install_runtime_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building EventNavigator_check_install_runtime.make"; \
	  $(cmtexe) -tag=$(tags) $(EventNavigator_check_install_runtime_extratags) build constituent_config -out=$(cmt_local_EventNavigator_check_install_runtime_makefile) EventNavigator_check_install_runtime; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_EventNavigator_check_install_runtime_makefile) : $(EventNavigator_check_install_runtime_dependencies) build_library_links
	$(echo) "(constituents.make) Building EventNavigator_check_install_runtime.make"; \
	  $(cmtexe) -f=$(bin)EventNavigator_check_install_runtime.in -tag=$(tags) $(EventNavigator_check_install_runtime_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_EventNavigator_check_install_runtime_makefile) EventNavigator_check_install_runtime
else
$(cmt_local_EventNavigator_check_install_runtime_makefile) : $(EventNavigator_check_install_runtime_dependencies) $(cmt_build_library_linksstamp) $(bin)EventNavigator_check_install_runtime.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_EventNavigator_check_install_runtime) ] || \
	  [ ! -f $(cmt_final_setup_EventNavigator_check_install_runtime) ] || \
	  $(not_EventNavigator_check_install_runtime_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building EventNavigator_check_install_runtime.make"; \
	  $(cmtexe) -f=$(bin)EventNavigator_check_install_runtime.in -tag=$(tags) $(EventNavigator_check_install_runtime_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_EventNavigator_check_install_runtime_makefile) EventNavigator_check_install_runtime; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(EventNavigator_check_install_runtime_extratags) build constituent_makefile -out=$(cmt_local_EventNavigator_check_install_runtime_makefile) EventNavigator_check_install_runtime

EventNavigator_check_install_runtime :: $(EventNavigator_check_install_runtime_dependencies) $(cmt_local_EventNavigator_check_install_runtime_makefile) dirs EventNavigator_check_install_runtimedirs
	$(echo) "(constituents.make) Starting EventNavigator_check_install_runtime"
	@if test -f $(cmt_local_EventNavigator_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_EventNavigator_check_install_runtime_makefile) EventNavigator_check_install_runtime; \
	  fi
#	@$(MAKE) -f $(cmt_local_EventNavigator_check_install_runtime_makefile) EventNavigator_check_install_runtime
	$(echo) "(constituents.make) EventNavigator_check_install_runtime done"

clean :: EventNavigator_check_install_runtimeclean

EventNavigator_check_install_runtimeclean :: $(EventNavigator_check_install_runtimeclean_dependencies) ##$(cmt_local_EventNavigator_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting EventNavigator_check_install_runtimeclean"
	@-if test -f $(cmt_local_EventNavigator_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_EventNavigator_check_install_runtime_makefile) EventNavigator_check_install_runtimeclean; \
	fi
	$(echo) "(constituents.make) EventNavigator_check_install_runtimeclean done"
#	@-$(MAKE) -f $(cmt_local_EventNavigator_check_install_runtime_makefile) EventNavigator_check_install_runtimeclean

##	  /bin/rm -f $(cmt_local_EventNavigator_check_install_runtime_makefile) $(bin)EventNavigator_check_install_runtime_dependencies.make

install :: EventNavigator_check_install_runtimeinstall

EventNavigator_check_install_runtimeinstall :: $(EventNavigator_check_install_runtime_dependencies) $(cmt_local_EventNavigator_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting install EventNavigator_check_install_runtime"
	@-$(MAKE) -f $(cmt_local_EventNavigator_check_install_runtime_makefile) install
	$(echo) "(constituents.make) install EventNavigator_check_install_runtime done"

uninstall : EventNavigator_check_install_runtimeuninstall

$(foreach d,$(EventNavigator_check_install_runtime_dependencies),$(eval $(d)uninstall_dependencies += EventNavigator_check_install_runtimeuninstall))

EventNavigator_check_install_runtimeuninstall : $(EventNavigator_check_install_runtimeuninstall_dependencies) ##$(cmt_local_EventNavigator_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting uninstall EventNavigator_check_install_runtime"
	@if test -f $(cmt_local_EventNavigator_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_EventNavigator_check_install_runtime_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_EventNavigator_check_install_runtime_makefile) uninstall
	$(echo) "(constituents.make) uninstall EventNavigator_check_install_runtime done"

remove_library_links :: EventNavigator_check_install_runtimeuninstall

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ EventNavigator_check_install_runtime"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ EventNavigator_check_install_runtime done"
endif

#-- end of constituent ------
#-- start of constituent ------

cmt_Gam4pikpAlg_check_install_runtime_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_Gam4pikpAlg_check_install_runtime_has_target_tag

#cmt_local_tagfile_Gam4pikpAlg_check_install_runtime = $(TestRelease_tag)_Gam4pikpAlg_check_install_runtime.make
cmt_local_tagfile_Gam4pikpAlg_check_install_runtime = $(bin)$(TestRelease_tag)_Gam4pikpAlg_check_install_runtime.make
cmt_local_setup_Gam4pikpAlg_check_install_runtime = $(bin)setup_Gam4pikpAlg_check_install_runtime$$$$.make
cmt_final_setup_Gam4pikpAlg_check_install_runtime = $(bin)setup_Gam4pikpAlg_check_install_runtime.make
#cmt_final_setup_Gam4pikpAlg_check_install_runtime = $(bin)TestRelease_Gam4pikpAlg_check_install_runtimesetup.make
cmt_local_Gam4pikpAlg_check_install_runtime_makefile = $(bin)Gam4pikpAlg_check_install_runtime.make

Gam4pikpAlg_check_install_runtime_extratags = -tag_add=target_Gam4pikpAlg_check_install_runtime

#$(cmt_local_tagfile_Gam4pikpAlg_check_install_runtime) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_Gam4pikpAlg_check_install_runtime) ::
else
$(cmt_local_tagfile_Gam4pikpAlg_check_install_runtime) :
endif
	$(echo) "(constituents.make) Rebuilding $@"; \
	  if test -f $(cmt_local_tagfile_Gam4pikpAlg_check_install_runtime); then /bin/rm -f $(cmt_local_tagfile_Gam4pikpAlg_check_install_runtime); fi ; \
	  $(cmtexe) -tag=$(tags) $(Gam4pikpAlg_check_install_runtime_extratags) build tag_makefile >>$(cmt_local_tagfile_Gam4pikpAlg_check_install_runtime)
	$(echo) "(constituents.make) Rebuilding $(cmt_final_setup_Gam4pikpAlg_check_install_runtime)"; \
	  test ! -f $(cmt_local_setup_Gam4pikpAlg_check_install_runtime) || \rm -f $(cmt_local_setup_Gam4pikpAlg_check_install_runtime); \
	  trap '\rm -f $(cmt_local_setup_Gam4pikpAlg_check_install_runtime)' 0 1 2 15; \
	  $(cmtexe) -tag=$(tags) $(Gam4pikpAlg_check_install_runtime_extratags) show setup >$(cmt_local_setup_Gam4pikpAlg_check_install_runtime) && \
	  if [ -f $(cmt_final_setup_Gam4pikpAlg_check_install_runtime) ] && \
	    \cmp -s $(cmt_final_setup_Gam4pikpAlg_check_install_runtime) $(cmt_local_setup_Gam4pikpAlg_check_install_runtime); then \
	    \rm $(cmt_local_setup_Gam4pikpAlg_check_install_runtime); else \
	    \mv -f $(cmt_local_setup_Gam4pikpAlg_check_install_runtime) $(cmt_final_setup_Gam4pikpAlg_check_install_runtime); fi

else

#cmt_local_tagfile_Gam4pikpAlg_check_install_runtime = $(TestRelease_tag).make
cmt_local_tagfile_Gam4pikpAlg_check_install_runtime = $(bin)$(TestRelease_tag).make
cmt_final_setup_Gam4pikpAlg_check_install_runtime = $(bin)setup.make
#cmt_final_setup_Gam4pikpAlg_check_install_runtime = $(bin)TestReleasesetup.make
cmt_local_Gam4pikpAlg_check_install_runtime_makefile = $(bin)Gam4pikpAlg_check_install_runtime.make

endif

not_Gam4pikpAlg_check_install_runtime_dependencies = { n=0; for p in $?; do m=0; for d in $(Gam4pikpAlg_check_install_runtime_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
Gam4pikpAlg_check_install_runtimedirs :
	@if test ! -d $(bin)Gam4pikpAlg_check_install_runtime; then $(mkdir) -p $(bin)Gam4pikpAlg_check_install_runtime; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)Gam4pikpAlg_check_install_runtime
else
Gam4pikpAlg_check_install_runtimedirs : ;
endif

#ifndef QUICK
#ifdef STRUCTURED_OUTPUT
# Gam4pikpAlg_check_install_runtimedirs ::
#	@if test ! -d $(bin)Gam4pikpAlg_check_install_runtime; then $(mkdir) -p $(bin)Gam4pikpAlg_check_install_runtime; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)Gam4pikpAlg_check_install_runtime
#
#$(cmt_local_Gam4pikpAlg_check_install_runtime_makefile) :: $(Gam4pikpAlg_check_install_runtime_dependencies) $(cmt_local_tagfile_Gam4pikpAlg_check_install_runtime) build_library_links dirs Gam4pikpAlg_check_install_runtimedirs
#else
#$(cmt_local_Gam4pikpAlg_check_install_runtime_makefile) :: $(Gam4pikpAlg_check_install_runtime_dependencies) $(cmt_local_tagfile_Gam4pikpAlg_check_install_runtime) build_library_links dirs
#endif
#else
#$(cmt_local_Gam4pikpAlg_check_install_runtime_makefile) :: $(cmt_local_tagfile_Gam4pikpAlg_check_install_runtime)
#endif

ifdef cmt_Gam4pikpAlg_check_install_runtime_has_target_tag

ifndef QUICK
$(cmt_local_Gam4pikpAlg_check_install_runtime_makefile) : $(Gam4pikpAlg_check_install_runtime_dependencies) build_library_links
	$(echo) "(constituents.make) Building Gam4pikpAlg_check_install_runtime.make"; \
	  $(cmtexe) -tag=$(tags) $(Gam4pikpAlg_check_install_runtime_extratags) build constituent_config -out=$(cmt_local_Gam4pikpAlg_check_install_runtime_makefile) Gam4pikpAlg_check_install_runtime
else
$(cmt_local_Gam4pikpAlg_check_install_runtime_makefile) : $(Gam4pikpAlg_check_install_runtime_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_Gam4pikpAlg_check_install_runtime) ] || \
	  [ ! -f $(cmt_final_setup_Gam4pikpAlg_check_install_runtime) ] || \
	  $(not_Gam4pikpAlg_check_install_runtime_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building Gam4pikpAlg_check_install_runtime.make"; \
	  $(cmtexe) -tag=$(tags) $(Gam4pikpAlg_check_install_runtime_extratags) build constituent_config -out=$(cmt_local_Gam4pikpAlg_check_install_runtime_makefile) Gam4pikpAlg_check_install_runtime; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_Gam4pikpAlg_check_install_runtime_makefile) : $(Gam4pikpAlg_check_install_runtime_dependencies) build_library_links
	$(echo) "(constituents.make) Building Gam4pikpAlg_check_install_runtime.make"; \
	  $(cmtexe) -f=$(bin)Gam4pikpAlg_check_install_runtime.in -tag=$(tags) $(Gam4pikpAlg_check_install_runtime_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_Gam4pikpAlg_check_install_runtime_makefile) Gam4pikpAlg_check_install_runtime
else
$(cmt_local_Gam4pikpAlg_check_install_runtime_makefile) : $(Gam4pikpAlg_check_install_runtime_dependencies) $(cmt_build_library_linksstamp) $(bin)Gam4pikpAlg_check_install_runtime.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_Gam4pikpAlg_check_install_runtime) ] || \
	  [ ! -f $(cmt_final_setup_Gam4pikpAlg_check_install_runtime) ] || \
	  $(not_Gam4pikpAlg_check_install_runtime_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building Gam4pikpAlg_check_install_runtime.make"; \
	  $(cmtexe) -f=$(bin)Gam4pikpAlg_check_install_runtime.in -tag=$(tags) $(Gam4pikpAlg_check_install_runtime_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_Gam4pikpAlg_check_install_runtime_makefile) Gam4pikpAlg_check_install_runtime; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(Gam4pikpAlg_check_install_runtime_extratags) build constituent_makefile -out=$(cmt_local_Gam4pikpAlg_check_install_runtime_makefile) Gam4pikpAlg_check_install_runtime

Gam4pikpAlg_check_install_runtime :: $(Gam4pikpAlg_check_install_runtime_dependencies) $(cmt_local_Gam4pikpAlg_check_install_runtime_makefile) dirs Gam4pikpAlg_check_install_runtimedirs
	$(echo) "(constituents.make) Starting Gam4pikpAlg_check_install_runtime"
	@if test -f $(cmt_local_Gam4pikpAlg_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_Gam4pikpAlg_check_install_runtime_makefile) Gam4pikpAlg_check_install_runtime; \
	  fi
#	@$(MAKE) -f $(cmt_local_Gam4pikpAlg_check_install_runtime_makefile) Gam4pikpAlg_check_install_runtime
	$(echo) "(constituents.make) Gam4pikpAlg_check_install_runtime done"

clean :: Gam4pikpAlg_check_install_runtimeclean

Gam4pikpAlg_check_install_runtimeclean :: $(Gam4pikpAlg_check_install_runtimeclean_dependencies) ##$(cmt_local_Gam4pikpAlg_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting Gam4pikpAlg_check_install_runtimeclean"
	@-if test -f $(cmt_local_Gam4pikpAlg_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_Gam4pikpAlg_check_install_runtime_makefile) Gam4pikpAlg_check_install_runtimeclean; \
	fi
	$(echo) "(constituents.make) Gam4pikpAlg_check_install_runtimeclean done"
#	@-$(MAKE) -f $(cmt_local_Gam4pikpAlg_check_install_runtime_makefile) Gam4pikpAlg_check_install_runtimeclean

##	  /bin/rm -f $(cmt_local_Gam4pikpAlg_check_install_runtime_makefile) $(bin)Gam4pikpAlg_check_install_runtime_dependencies.make

install :: Gam4pikpAlg_check_install_runtimeinstall

Gam4pikpAlg_check_install_runtimeinstall :: $(Gam4pikpAlg_check_install_runtime_dependencies) $(cmt_local_Gam4pikpAlg_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting install Gam4pikpAlg_check_install_runtime"
	@-$(MAKE) -f $(cmt_local_Gam4pikpAlg_check_install_runtime_makefile) install
	$(echo) "(constituents.make) install Gam4pikpAlg_check_install_runtime done"

uninstall : Gam4pikpAlg_check_install_runtimeuninstall

$(foreach d,$(Gam4pikpAlg_check_install_runtime_dependencies),$(eval $(d)uninstall_dependencies += Gam4pikpAlg_check_install_runtimeuninstall))

Gam4pikpAlg_check_install_runtimeuninstall : $(Gam4pikpAlg_check_install_runtimeuninstall_dependencies) ##$(cmt_local_Gam4pikpAlg_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting uninstall Gam4pikpAlg_check_install_runtime"
	@if test -f $(cmt_local_Gam4pikpAlg_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_Gam4pikpAlg_check_install_runtime_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_Gam4pikpAlg_check_install_runtime_makefile) uninstall
	$(echo) "(constituents.make) uninstall Gam4pikpAlg_check_install_runtime done"

remove_library_links :: Gam4pikpAlg_check_install_runtimeuninstall

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ Gam4pikpAlg_check_install_runtime"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ Gam4pikpAlg_check_install_runtime done"
endif

#-- end of constituent ------
#-- start of constituent ------

cmt_BesEvtGen_check_install_runtime_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_BesEvtGen_check_install_runtime_has_target_tag

#cmt_local_tagfile_BesEvtGen_check_install_runtime = $(TestRelease_tag)_BesEvtGen_check_install_runtime.make
cmt_local_tagfile_BesEvtGen_check_install_runtime = $(bin)$(TestRelease_tag)_BesEvtGen_check_install_runtime.make
cmt_local_setup_BesEvtGen_check_install_runtime = $(bin)setup_BesEvtGen_check_install_runtime$$$$.make
cmt_final_setup_BesEvtGen_check_install_runtime = $(bin)setup_BesEvtGen_check_install_runtime.make
#cmt_final_setup_BesEvtGen_check_install_runtime = $(bin)TestRelease_BesEvtGen_check_install_runtimesetup.make
cmt_local_BesEvtGen_check_install_runtime_makefile = $(bin)BesEvtGen_check_install_runtime.make

BesEvtGen_check_install_runtime_extratags = -tag_add=target_BesEvtGen_check_install_runtime

#$(cmt_local_tagfile_BesEvtGen_check_install_runtime) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_BesEvtGen_check_install_runtime) ::
else
$(cmt_local_tagfile_BesEvtGen_check_install_runtime) :
endif
	$(echo) "(constituents.make) Rebuilding $@"; \
	  if test -f $(cmt_local_tagfile_BesEvtGen_check_install_runtime); then /bin/rm -f $(cmt_local_tagfile_BesEvtGen_check_install_runtime); fi ; \
	  $(cmtexe) -tag=$(tags) $(BesEvtGen_check_install_runtime_extratags) build tag_makefile >>$(cmt_local_tagfile_BesEvtGen_check_install_runtime)
	$(echo) "(constituents.make) Rebuilding $(cmt_final_setup_BesEvtGen_check_install_runtime)"; \
	  test ! -f $(cmt_local_setup_BesEvtGen_check_install_runtime) || \rm -f $(cmt_local_setup_BesEvtGen_check_install_runtime); \
	  trap '\rm -f $(cmt_local_setup_BesEvtGen_check_install_runtime)' 0 1 2 15; \
	  $(cmtexe) -tag=$(tags) $(BesEvtGen_check_install_runtime_extratags) show setup >$(cmt_local_setup_BesEvtGen_check_install_runtime) && \
	  if [ -f $(cmt_final_setup_BesEvtGen_check_install_runtime) ] && \
	    \cmp -s $(cmt_final_setup_BesEvtGen_check_install_runtime) $(cmt_local_setup_BesEvtGen_check_install_runtime); then \
	    \rm $(cmt_local_setup_BesEvtGen_check_install_runtime); else \
	    \mv -f $(cmt_local_setup_BesEvtGen_check_install_runtime) $(cmt_final_setup_BesEvtGen_check_install_runtime); fi

else

#cmt_local_tagfile_BesEvtGen_check_install_runtime = $(TestRelease_tag).make
cmt_local_tagfile_BesEvtGen_check_install_runtime = $(bin)$(TestRelease_tag).make
cmt_final_setup_BesEvtGen_check_install_runtime = $(bin)setup.make
#cmt_final_setup_BesEvtGen_check_install_runtime = $(bin)TestReleasesetup.make
cmt_local_BesEvtGen_check_install_runtime_makefile = $(bin)BesEvtGen_check_install_runtime.make

endif

not_BesEvtGen_check_install_runtime_dependencies = { n=0; for p in $?; do m=0; for d in $(BesEvtGen_check_install_runtime_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
BesEvtGen_check_install_runtimedirs :
	@if test ! -d $(bin)BesEvtGen_check_install_runtime; then $(mkdir) -p $(bin)BesEvtGen_check_install_runtime; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)BesEvtGen_check_install_runtime
else
BesEvtGen_check_install_runtimedirs : ;
endif

#ifndef QUICK
#ifdef STRUCTURED_OUTPUT
# BesEvtGen_check_install_runtimedirs ::
#	@if test ! -d $(bin)BesEvtGen_check_install_runtime; then $(mkdir) -p $(bin)BesEvtGen_check_install_runtime; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)BesEvtGen_check_install_runtime
#
#$(cmt_local_BesEvtGen_check_install_runtime_makefile) :: $(BesEvtGen_check_install_runtime_dependencies) $(cmt_local_tagfile_BesEvtGen_check_install_runtime) build_library_links dirs BesEvtGen_check_install_runtimedirs
#else
#$(cmt_local_BesEvtGen_check_install_runtime_makefile) :: $(BesEvtGen_check_install_runtime_dependencies) $(cmt_local_tagfile_BesEvtGen_check_install_runtime) build_library_links dirs
#endif
#else
#$(cmt_local_BesEvtGen_check_install_runtime_makefile) :: $(cmt_local_tagfile_BesEvtGen_check_install_runtime)
#endif

ifdef cmt_BesEvtGen_check_install_runtime_has_target_tag

ifndef QUICK
$(cmt_local_BesEvtGen_check_install_runtime_makefile) : $(BesEvtGen_check_install_runtime_dependencies) build_library_links
	$(echo) "(constituents.make) Building BesEvtGen_check_install_runtime.make"; \
	  $(cmtexe) -tag=$(tags) $(BesEvtGen_check_install_runtime_extratags) build constituent_config -out=$(cmt_local_BesEvtGen_check_install_runtime_makefile) BesEvtGen_check_install_runtime
else
$(cmt_local_BesEvtGen_check_install_runtime_makefile) : $(BesEvtGen_check_install_runtime_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_BesEvtGen_check_install_runtime) ] || \
	  [ ! -f $(cmt_final_setup_BesEvtGen_check_install_runtime) ] || \
	  $(not_BesEvtGen_check_install_runtime_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building BesEvtGen_check_install_runtime.make"; \
	  $(cmtexe) -tag=$(tags) $(BesEvtGen_check_install_runtime_extratags) build constituent_config -out=$(cmt_local_BesEvtGen_check_install_runtime_makefile) BesEvtGen_check_install_runtime; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_BesEvtGen_check_install_runtime_makefile) : $(BesEvtGen_check_install_runtime_dependencies) build_library_links
	$(echo) "(constituents.make) Building BesEvtGen_check_install_runtime.make"; \
	  $(cmtexe) -f=$(bin)BesEvtGen_check_install_runtime.in -tag=$(tags) $(BesEvtGen_check_install_runtime_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_BesEvtGen_check_install_runtime_makefile) BesEvtGen_check_install_runtime
else
$(cmt_local_BesEvtGen_check_install_runtime_makefile) : $(BesEvtGen_check_install_runtime_dependencies) $(cmt_build_library_linksstamp) $(bin)BesEvtGen_check_install_runtime.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_BesEvtGen_check_install_runtime) ] || \
	  [ ! -f $(cmt_final_setup_BesEvtGen_check_install_runtime) ] || \
	  $(not_BesEvtGen_check_install_runtime_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building BesEvtGen_check_install_runtime.make"; \
	  $(cmtexe) -f=$(bin)BesEvtGen_check_install_runtime.in -tag=$(tags) $(BesEvtGen_check_install_runtime_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_BesEvtGen_check_install_runtime_makefile) BesEvtGen_check_install_runtime; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(BesEvtGen_check_install_runtime_extratags) build constituent_makefile -out=$(cmt_local_BesEvtGen_check_install_runtime_makefile) BesEvtGen_check_install_runtime

BesEvtGen_check_install_runtime :: $(BesEvtGen_check_install_runtime_dependencies) $(cmt_local_BesEvtGen_check_install_runtime_makefile) dirs BesEvtGen_check_install_runtimedirs
	$(echo) "(constituents.make) Starting BesEvtGen_check_install_runtime"
	@if test -f $(cmt_local_BesEvtGen_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_BesEvtGen_check_install_runtime_makefile) BesEvtGen_check_install_runtime; \
	  fi
#	@$(MAKE) -f $(cmt_local_BesEvtGen_check_install_runtime_makefile) BesEvtGen_check_install_runtime
	$(echo) "(constituents.make) BesEvtGen_check_install_runtime done"

clean :: BesEvtGen_check_install_runtimeclean

BesEvtGen_check_install_runtimeclean :: $(BesEvtGen_check_install_runtimeclean_dependencies) ##$(cmt_local_BesEvtGen_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting BesEvtGen_check_install_runtimeclean"
	@-if test -f $(cmt_local_BesEvtGen_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_BesEvtGen_check_install_runtime_makefile) BesEvtGen_check_install_runtimeclean; \
	fi
	$(echo) "(constituents.make) BesEvtGen_check_install_runtimeclean done"
#	@-$(MAKE) -f $(cmt_local_BesEvtGen_check_install_runtime_makefile) BesEvtGen_check_install_runtimeclean

##	  /bin/rm -f $(cmt_local_BesEvtGen_check_install_runtime_makefile) $(bin)BesEvtGen_check_install_runtime_dependencies.make

install :: BesEvtGen_check_install_runtimeinstall

BesEvtGen_check_install_runtimeinstall :: $(BesEvtGen_check_install_runtime_dependencies) $(cmt_local_BesEvtGen_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting install BesEvtGen_check_install_runtime"
	@-$(MAKE) -f $(cmt_local_BesEvtGen_check_install_runtime_makefile) install
	$(echo) "(constituents.make) install BesEvtGen_check_install_runtime done"

uninstall : BesEvtGen_check_install_runtimeuninstall

$(foreach d,$(BesEvtGen_check_install_runtime_dependencies),$(eval $(d)uninstall_dependencies += BesEvtGen_check_install_runtimeuninstall))

BesEvtGen_check_install_runtimeuninstall : $(BesEvtGen_check_install_runtimeuninstall_dependencies) ##$(cmt_local_BesEvtGen_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting uninstall BesEvtGen_check_install_runtime"
	@if test -f $(cmt_local_BesEvtGen_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_BesEvtGen_check_install_runtime_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_BesEvtGen_check_install_runtime_makefile) uninstall
	$(echo) "(constituents.make) uninstall BesEvtGen_check_install_runtime done"

remove_library_links :: BesEvtGen_check_install_runtimeuninstall

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ BesEvtGen_check_install_runtime"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ BesEvtGen_check_install_runtime done"
endif

#-- end of constituent ------
#-- start of constituent ------

cmt_Bhwide_check_install_runtime_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_Bhwide_check_install_runtime_has_target_tag

#cmt_local_tagfile_Bhwide_check_install_runtime = $(TestRelease_tag)_Bhwide_check_install_runtime.make
cmt_local_tagfile_Bhwide_check_install_runtime = $(bin)$(TestRelease_tag)_Bhwide_check_install_runtime.make
cmt_local_setup_Bhwide_check_install_runtime = $(bin)setup_Bhwide_check_install_runtime$$$$.make
cmt_final_setup_Bhwide_check_install_runtime = $(bin)setup_Bhwide_check_install_runtime.make
#cmt_final_setup_Bhwide_check_install_runtime = $(bin)TestRelease_Bhwide_check_install_runtimesetup.make
cmt_local_Bhwide_check_install_runtime_makefile = $(bin)Bhwide_check_install_runtime.make

Bhwide_check_install_runtime_extratags = -tag_add=target_Bhwide_check_install_runtime

#$(cmt_local_tagfile_Bhwide_check_install_runtime) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_Bhwide_check_install_runtime) ::
else
$(cmt_local_tagfile_Bhwide_check_install_runtime) :
endif
	$(echo) "(constituents.make) Rebuilding $@"; \
	  if test -f $(cmt_local_tagfile_Bhwide_check_install_runtime); then /bin/rm -f $(cmt_local_tagfile_Bhwide_check_install_runtime); fi ; \
	  $(cmtexe) -tag=$(tags) $(Bhwide_check_install_runtime_extratags) build tag_makefile >>$(cmt_local_tagfile_Bhwide_check_install_runtime)
	$(echo) "(constituents.make) Rebuilding $(cmt_final_setup_Bhwide_check_install_runtime)"; \
	  test ! -f $(cmt_local_setup_Bhwide_check_install_runtime) || \rm -f $(cmt_local_setup_Bhwide_check_install_runtime); \
	  trap '\rm -f $(cmt_local_setup_Bhwide_check_install_runtime)' 0 1 2 15; \
	  $(cmtexe) -tag=$(tags) $(Bhwide_check_install_runtime_extratags) show setup >$(cmt_local_setup_Bhwide_check_install_runtime) && \
	  if [ -f $(cmt_final_setup_Bhwide_check_install_runtime) ] && \
	    \cmp -s $(cmt_final_setup_Bhwide_check_install_runtime) $(cmt_local_setup_Bhwide_check_install_runtime); then \
	    \rm $(cmt_local_setup_Bhwide_check_install_runtime); else \
	    \mv -f $(cmt_local_setup_Bhwide_check_install_runtime) $(cmt_final_setup_Bhwide_check_install_runtime); fi

else

#cmt_local_tagfile_Bhwide_check_install_runtime = $(TestRelease_tag).make
cmt_local_tagfile_Bhwide_check_install_runtime = $(bin)$(TestRelease_tag).make
cmt_final_setup_Bhwide_check_install_runtime = $(bin)setup.make
#cmt_final_setup_Bhwide_check_install_runtime = $(bin)TestReleasesetup.make
cmt_local_Bhwide_check_install_runtime_makefile = $(bin)Bhwide_check_install_runtime.make

endif

not_Bhwide_check_install_runtime_dependencies = { n=0; for p in $?; do m=0; for d in $(Bhwide_check_install_runtime_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
Bhwide_check_install_runtimedirs :
	@if test ! -d $(bin)Bhwide_check_install_runtime; then $(mkdir) -p $(bin)Bhwide_check_install_runtime; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)Bhwide_check_install_runtime
else
Bhwide_check_install_runtimedirs : ;
endif

#ifndef QUICK
#ifdef STRUCTURED_OUTPUT
# Bhwide_check_install_runtimedirs ::
#	@if test ! -d $(bin)Bhwide_check_install_runtime; then $(mkdir) -p $(bin)Bhwide_check_install_runtime; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)Bhwide_check_install_runtime
#
#$(cmt_local_Bhwide_check_install_runtime_makefile) :: $(Bhwide_check_install_runtime_dependencies) $(cmt_local_tagfile_Bhwide_check_install_runtime) build_library_links dirs Bhwide_check_install_runtimedirs
#else
#$(cmt_local_Bhwide_check_install_runtime_makefile) :: $(Bhwide_check_install_runtime_dependencies) $(cmt_local_tagfile_Bhwide_check_install_runtime) build_library_links dirs
#endif
#else
#$(cmt_local_Bhwide_check_install_runtime_makefile) :: $(cmt_local_tagfile_Bhwide_check_install_runtime)
#endif

ifdef cmt_Bhwide_check_install_runtime_has_target_tag

ifndef QUICK
$(cmt_local_Bhwide_check_install_runtime_makefile) : $(Bhwide_check_install_runtime_dependencies) build_library_links
	$(echo) "(constituents.make) Building Bhwide_check_install_runtime.make"; \
	  $(cmtexe) -tag=$(tags) $(Bhwide_check_install_runtime_extratags) build constituent_config -out=$(cmt_local_Bhwide_check_install_runtime_makefile) Bhwide_check_install_runtime
else
$(cmt_local_Bhwide_check_install_runtime_makefile) : $(Bhwide_check_install_runtime_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_Bhwide_check_install_runtime) ] || \
	  [ ! -f $(cmt_final_setup_Bhwide_check_install_runtime) ] || \
	  $(not_Bhwide_check_install_runtime_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building Bhwide_check_install_runtime.make"; \
	  $(cmtexe) -tag=$(tags) $(Bhwide_check_install_runtime_extratags) build constituent_config -out=$(cmt_local_Bhwide_check_install_runtime_makefile) Bhwide_check_install_runtime; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_Bhwide_check_install_runtime_makefile) : $(Bhwide_check_install_runtime_dependencies) build_library_links
	$(echo) "(constituents.make) Building Bhwide_check_install_runtime.make"; \
	  $(cmtexe) -f=$(bin)Bhwide_check_install_runtime.in -tag=$(tags) $(Bhwide_check_install_runtime_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_Bhwide_check_install_runtime_makefile) Bhwide_check_install_runtime
else
$(cmt_local_Bhwide_check_install_runtime_makefile) : $(Bhwide_check_install_runtime_dependencies) $(cmt_build_library_linksstamp) $(bin)Bhwide_check_install_runtime.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_Bhwide_check_install_runtime) ] || \
	  [ ! -f $(cmt_final_setup_Bhwide_check_install_runtime) ] || \
	  $(not_Bhwide_check_install_runtime_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building Bhwide_check_install_runtime.make"; \
	  $(cmtexe) -f=$(bin)Bhwide_check_install_runtime.in -tag=$(tags) $(Bhwide_check_install_runtime_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_Bhwide_check_install_runtime_makefile) Bhwide_check_install_runtime; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(Bhwide_check_install_runtime_extratags) build constituent_makefile -out=$(cmt_local_Bhwide_check_install_runtime_makefile) Bhwide_check_install_runtime

Bhwide_check_install_runtime :: $(Bhwide_check_install_runtime_dependencies) $(cmt_local_Bhwide_check_install_runtime_makefile) dirs Bhwide_check_install_runtimedirs
	$(echo) "(constituents.make) Starting Bhwide_check_install_runtime"
	@if test -f $(cmt_local_Bhwide_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_Bhwide_check_install_runtime_makefile) Bhwide_check_install_runtime; \
	  fi
#	@$(MAKE) -f $(cmt_local_Bhwide_check_install_runtime_makefile) Bhwide_check_install_runtime
	$(echo) "(constituents.make) Bhwide_check_install_runtime done"

clean :: Bhwide_check_install_runtimeclean

Bhwide_check_install_runtimeclean :: $(Bhwide_check_install_runtimeclean_dependencies) ##$(cmt_local_Bhwide_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting Bhwide_check_install_runtimeclean"
	@-if test -f $(cmt_local_Bhwide_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_Bhwide_check_install_runtime_makefile) Bhwide_check_install_runtimeclean; \
	fi
	$(echo) "(constituents.make) Bhwide_check_install_runtimeclean done"
#	@-$(MAKE) -f $(cmt_local_Bhwide_check_install_runtime_makefile) Bhwide_check_install_runtimeclean

##	  /bin/rm -f $(cmt_local_Bhwide_check_install_runtime_makefile) $(bin)Bhwide_check_install_runtime_dependencies.make

install :: Bhwide_check_install_runtimeinstall

Bhwide_check_install_runtimeinstall :: $(Bhwide_check_install_runtime_dependencies) $(cmt_local_Bhwide_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting install Bhwide_check_install_runtime"
	@-$(MAKE) -f $(cmt_local_Bhwide_check_install_runtime_makefile) install
	$(echo) "(constituents.make) install Bhwide_check_install_runtime done"

uninstall : Bhwide_check_install_runtimeuninstall

$(foreach d,$(Bhwide_check_install_runtime_dependencies),$(eval $(d)uninstall_dependencies += Bhwide_check_install_runtimeuninstall))

Bhwide_check_install_runtimeuninstall : $(Bhwide_check_install_runtimeuninstall_dependencies) ##$(cmt_local_Bhwide_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting uninstall Bhwide_check_install_runtime"
	@if test -f $(cmt_local_Bhwide_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_Bhwide_check_install_runtime_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_Bhwide_check_install_runtime_makefile) uninstall
	$(echo) "(constituents.make) uninstall Bhwide_check_install_runtime done"

remove_library_links :: Bhwide_check_install_runtimeuninstall

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ Bhwide_check_install_runtime"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ Bhwide_check_install_runtime done"
endif

#-- end of constituent ------
#-- start of constituent ------

cmt_Babayaga_check_install_runtime_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_Babayaga_check_install_runtime_has_target_tag

#cmt_local_tagfile_Babayaga_check_install_runtime = $(TestRelease_tag)_Babayaga_check_install_runtime.make
cmt_local_tagfile_Babayaga_check_install_runtime = $(bin)$(TestRelease_tag)_Babayaga_check_install_runtime.make
cmt_local_setup_Babayaga_check_install_runtime = $(bin)setup_Babayaga_check_install_runtime$$$$.make
cmt_final_setup_Babayaga_check_install_runtime = $(bin)setup_Babayaga_check_install_runtime.make
#cmt_final_setup_Babayaga_check_install_runtime = $(bin)TestRelease_Babayaga_check_install_runtimesetup.make
cmt_local_Babayaga_check_install_runtime_makefile = $(bin)Babayaga_check_install_runtime.make

Babayaga_check_install_runtime_extratags = -tag_add=target_Babayaga_check_install_runtime

#$(cmt_local_tagfile_Babayaga_check_install_runtime) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_Babayaga_check_install_runtime) ::
else
$(cmt_local_tagfile_Babayaga_check_install_runtime) :
endif
	$(echo) "(constituents.make) Rebuilding $@"; \
	  if test -f $(cmt_local_tagfile_Babayaga_check_install_runtime); then /bin/rm -f $(cmt_local_tagfile_Babayaga_check_install_runtime); fi ; \
	  $(cmtexe) -tag=$(tags) $(Babayaga_check_install_runtime_extratags) build tag_makefile >>$(cmt_local_tagfile_Babayaga_check_install_runtime)
	$(echo) "(constituents.make) Rebuilding $(cmt_final_setup_Babayaga_check_install_runtime)"; \
	  test ! -f $(cmt_local_setup_Babayaga_check_install_runtime) || \rm -f $(cmt_local_setup_Babayaga_check_install_runtime); \
	  trap '\rm -f $(cmt_local_setup_Babayaga_check_install_runtime)' 0 1 2 15; \
	  $(cmtexe) -tag=$(tags) $(Babayaga_check_install_runtime_extratags) show setup >$(cmt_local_setup_Babayaga_check_install_runtime) && \
	  if [ -f $(cmt_final_setup_Babayaga_check_install_runtime) ] && \
	    \cmp -s $(cmt_final_setup_Babayaga_check_install_runtime) $(cmt_local_setup_Babayaga_check_install_runtime); then \
	    \rm $(cmt_local_setup_Babayaga_check_install_runtime); else \
	    \mv -f $(cmt_local_setup_Babayaga_check_install_runtime) $(cmt_final_setup_Babayaga_check_install_runtime); fi

else

#cmt_local_tagfile_Babayaga_check_install_runtime = $(TestRelease_tag).make
cmt_local_tagfile_Babayaga_check_install_runtime = $(bin)$(TestRelease_tag).make
cmt_final_setup_Babayaga_check_install_runtime = $(bin)setup.make
#cmt_final_setup_Babayaga_check_install_runtime = $(bin)TestReleasesetup.make
cmt_local_Babayaga_check_install_runtime_makefile = $(bin)Babayaga_check_install_runtime.make

endif

not_Babayaga_check_install_runtime_dependencies = { n=0; for p in $?; do m=0; for d in $(Babayaga_check_install_runtime_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
Babayaga_check_install_runtimedirs :
	@if test ! -d $(bin)Babayaga_check_install_runtime; then $(mkdir) -p $(bin)Babayaga_check_install_runtime; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)Babayaga_check_install_runtime
else
Babayaga_check_install_runtimedirs : ;
endif

#ifndef QUICK
#ifdef STRUCTURED_OUTPUT
# Babayaga_check_install_runtimedirs ::
#	@if test ! -d $(bin)Babayaga_check_install_runtime; then $(mkdir) -p $(bin)Babayaga_check_install_runtime; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)Babayaga_check_install_runtime
#
#$(cmt_local_Babayaga_check_install_runtime_makefile) :: $(Babayaga_check_install_runtime_dependencies) $(cmt_local_tagfile_Babayaga_check_install_runtime) build_library_links dirs Babayaga_check_install_runtimedirs
#else
#$(cmt_local_Babayaga_check_install_runtime_makefile) :: $(Babayaga_check_install_runtime_dependencies) $(cmt_local_tagfile_Babayaga_check_install_runtime) build_library_links dirs
#endif
#else
#$(cmt_local_Babayaga_check_install_runtime_makefile) :: $(cmt_local_tagfile_Babayaga_check_install_runtime)
#endif

ifdef cmt_Babayaga_check_install_runtime_has_target_tag

ifndef QUICK
$(cmt_local_Babayaga_check_install_runtime_makefile) : $(Babayaga_check_install_runtime_dependencies) build_library_links
	$(echo) "(constituents.make) Building Babayaga_check_install_runtime.make"; \
	  $(cmtexe) -tag=$(tags) $(Babayaga_check_install_runtime_extratags) build constituent_config -out=$(cmt_local_Babayaga_check_install_runtime_makefile) Babayaga_check_install_runtime
else
$(cmt_local_Babayaga_check_install_runtime_makefile) : $(Babayaga_check_install_runtime_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_Babayaga_check_install_runtime) ] || \
	  [ ! -f $(cmt_final_setup_Babayaga_check_install_runtime) ] || \
	  $(not_Babayaga_check_install_runtime_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building Babayaga_check_install_runtime.make"; \
	  $(cmtexe) -tag=$(tags) $(Babayaga_check_install_runtime_extratags) build constituent_config -out=$(cmt_local_Babayaga_check_install_runtime_makefile) Babayaga_check_install_runtime; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_Babayaga_check_install_runtime_makefile) : $(Babayaga_check_install_runtime_dependencies) build_library_links
	$(echo) "(constituents.make) Building Babayaga_check_install_runtime.make"; \
	  $(cmtexe) -f=$(bin)Babayaga_check_install_runtime.in -tag=$(tags) $(Babayaga_check_install_runtime_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_Babayaga_check_install_runtime_makefile) Babayaga_check_install_runtime
else
$(cmt_local_Babayaga_check_install_runtime_makefile) : $(Babayaga_check_install_runtime_dependencies) $(cmt_build_library_linksstamp) $(bin)Babayaga_check_install_runtime.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_Babayaga_check_install_runtime) ] || \
	  [ ! -f $(cmt_final_setup_Babayaga_check_install_runtime) ] || \
	  $(not_Babayaga_check_install_runtime_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building Babayaga_check_install_runtime.make"; \
	  $(cmtexe) -f=$(bin)Babayaga_check_install_runtime.in -tag=$(tags) $(Babayaga_check_install_runtime_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_Babayaga_check_install_runtime_makefile) Babayaga_check_install_runtime; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(Babayaga_check_install_runtime_extratags) build constituent_makefile -out=$(cmt_local_Babayaga_check_install_runtime_makefile) Babayaga_check_install_runtime

Babayaga_check_install_runtime :: $(Babayaga_check_install_runtime_dependencies) $(cmt_local_Babayaga_check_install_runtime_makefile) dirs Babayaga_check_install_runtimedirs
	$(echo) "(constituents.make) Starting Babayaga_check_install_runtime"
	@if test -f $(cmt_local_Babayaga_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_Babayaga_check_install_runtime_makefile) Babayaga_check_install_runtime; \
	  fi
#	@$(MAKE) -f $(cmt_local_Babayaga_check_install_runtime_makefile) Babayaga_check_install_runtime
	$(echo) "(constituents.make) Babayaga_check_install_runtime done"

clean :: Babayaga_check_install_runtimeclean

Babayaga_check_install_runtimeclean :: $(Babayaga_check_install_runtimeclean_dependencies) ##$(cmt_local_Babayaga_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting Babayaga_check_install_runtimeclean"
	@-if test -f $(cmt_local_Babayaga_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_Babayaga_check_install_runtime_makefile) Babayaga_check_install_runtimeclean; \
	fi
	$(echo) "(constituents.make) Babayaga_check_install_runtimeclean done"
#	@-$(MAKE) -f $(cmt_local_Babayaga_check_install_runtime_makefile) Babayaga_check_install_runtimeclean

##	  /bin/rm -f $(cmt_local_Babayaga_check_install_runtime_makefile) $(bin)Babayaga_check_install_runtime_dependencies.make

install :: Babayaga_check_install_runtimeinstall

Babayaga_check_install_runtimeinstall :: $(Babayaga_check_install_runtime_dependencies) $(cmt_local_Babayaga_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting install Babayaga_check_install_runtime"
	@-$(MAKE) -f $(cmt_local_Babayaga_check_install_runtime_makefile) install
	$(echo) "(constituents.make) install Babayaga_check_install_runtime done"

uninstall : Babayaga_check_install_runtimeuninstall

$(foreach d,$(Babayaga_check_install_runtime_dependencies),$(eval $(d)uninstall_dependencies += Babayaga_check_install_runtimeuninstall))

Babayaga_check_install_runtimeuninstall : $(Babayaga_check_install_runtimeuninstall_dependencies) ##$(cmt_local_Babayaga_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting uninstall Babayaga_check_install_runtime"
	@if test -f $(cmt_local_Babayaga_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_Babayaga_check_install_runtime_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_Babayaga_check_install_runtime_makefile) uninstall
	$(echo) "(constituents.make) uninstall Babayaga_check_install_runtime done"

remove_library_links :: Babayaga_check_install_runtimeuninstall

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ Babayaga_check_install_runtime"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ Babayaga_check_install_runtime done"
endif

#-- end of constituent ------
#-- start of constituent ------

cmt_RunEventNumberAlg_check_install_runtime_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_RunEventNumberAlg_check_install_runtime_has_target_tag

#cmt_local_tagfile_RunEventNumberAlg_check_install_runtime = $(TestRelease_tag)_RunEventNumberAlg_check_install_runtime.make
cmt_local_tagfile_RunEventNumberAlg_check_install_runtime = $(bin)$(TestRelease_tag)_RunEventNumberAlg_check_install_runtime.make
cmt_local_setup_RunEventNumberAlg_check_install_runtime = $(bin)setup_RunEventNumberAlg_check_install_runtime$$$$.make
cmt_final_setup_RunEventNumberAlg_check_install_runtime = $(bin)setup_RunEventNumberAlg_check_install_runtime.make
#cmt_final_setup_RunEventNumberAlg_check_install_runtime = $(bin)TestRelease_RunEventNumberAlg_check_install_runtimesetup.make
cmt_local_RunEventNumberAlg_check_install_runtime_makefile = $(bin)RunEventNumberAlg_check_install_runtime.make

RunEventNumberAlg_check_install_runtime_extratags = -tag_add=target_RunEventNumberAlg_check_install_runtime

#$(cmt_local_tagfile_RunEventNumberAlg_check_install_runtime) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_RunEventNumberAlg_check_install_runtime) ::
else
$(cmt_local_tagfile_RunEventNumberAlg_check_install_runtime) :
endif
	$(echo) "(constituents.make) Rebuilding $@"; \
	  if test -f $(cmt_local_tagfile_RunEventNumberAlg_check_install_runtime); then /bin/rm -f $(cmt_local_tagfile_RunEventNumberAlg_check_install_runtime); fi ; \
	  $(cmtexe) -tag=$(tags) $(RunEventNumberAlg_check_install_runtime_extratags) build tag_makefile >>$(cmt_local_tagfile_RunEventNumberAlg_check_install_runtime)
	$(echo) "(constituents.make) Rebuilding $(cmt_final_setup_RunEventNumberAlg_check_install_runtime)"; \
	  test ! -f $(cmt_local_setup_RunEventNumberAlg_check_install_runtime) || \rm -f $(cmt_local_setup_RunEventNumberAlg_check_install_runtime); \
	  trap '\rm -f $(cmt_local_setup_RunEventNumberAlg_check_install_runtime)' 0 1 2 15; \
	  $(cmtexe) -tag=$(tags) $(RunEventNumberAlg_check_install_runtime_extratags) show setup >$(cmt_local_setup_RunEventNumberAlg_check_install_runtime) && \
	  if [ -f $(cmt_final_setup_RunEventNumberAlg_check_install_runtime) ] && \
	    \cmp -s $(cmt_final_setup_RunEventNumberAlg_check_install_runtime) $(cmt_local_setup_RunEventNumberAlg_check_install_runtime); then \
	    \rm $(cmt_local_setup_RunEventNumberAlg_check_install_runtime); else \
	    \mv -f $(cmt_local_setup_RunEventNumberAlg_check_install_runtime) $(cmt_final_setup_RunEventNumberAlg_check_install_runtime); fi

else

#cmt_local_tagfile_RunEventNumberAlg_check_install_runtime = $(TestRelease_tag).make
cmt_local_tagfile_RunEventNumberAlg_check_install_runtime = $(bin)$(TestRelease_tag).make
cmt_final_setup_RunEventNumberAlg_check_install_runtime = $(bin)setup.make
#cmt_final_setup_RunEventNumberAlg_check_install_runtime = $(bin)TestReleasesetup.make
cmt_local_RunEventNumberAlg_check_install_runtime_makefile = $(bin)RunEventNumberAlg_check_install_runtime.make

endif

not_RunEventNumberAlg_check_install_runtime_dependencies = { n=0; for p in $?; do m=0; for d in $(RunEventNumberAlg_check_install_runtime_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
RunEventNumberAlg_check_install_runtimedirs :
	@if test ! -d $(bin)RunEventNumberAlg_check_install_runtime; then $(mkdir) -p $(bin)RunEventNumberAlg_check_install_runtime; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)RunEventNumberAlg_check_install_runtime
else
RunEventNumberAlg_check_install_runtimedirs : ;
endif

#ifndef QUICK
#ifdef STRUCTURED_OUTPUT
# RunEventNumberAlg_check_install_runtimedirs ::
#	@if test ! -d $(bin)RunEventNumberAlg_check_install_runtime; then $(mkdir) -p $(bin)RunEventNumberAlg_check_install_runtime; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)RunEventNumberAlg_check_install_runtime
#
#$(cmt_local_RunEventNumberAlg_check_install_runtime_makefile) :: $(RunEventNumberAlg_check_install_runtime_dependencies) $(cmt_local_tagfile_RunEventNumberAlg_check_install_runtime) build_library_links dirs RunEventNumberAlg_check_install_runtimedirs
#else
#$(cmt_local_RunEventNumberAlg_check_install_runtime_makefile) :: $(RunEventNumberAlg_check_install_runtime_dependencies) $(cmt_local_tagfile_RunEventNumberAlg_check_install_runtime) build_library_links dirs
#endif
#else
#$(cmt_local_RunEventNumberAlg_check_install_runtime_makefile) :: $(cmt_local_tagfile_RunEventNumberAlg_check_install_runtime)
#endif

ifdef cmt_RunEventNumberAlg_check_install_runtime_has_target_tag

ifndef QUICK
$(cmt_local_RunEventNumberAlg_check_install_runtime_makefile) : $(RunEventNumberAlg_check_install_runtime_dependencies) build_library_links
	$(echo) "(constituents.make) Building RunEventNumberAlg_check_install_runtime.make"; \
	  $(cmtexe) -tag=$(tags) $(RunEventNumberAlg_check_install_runtime_extratags) build constituent_config -out=$(cmt_local_RunEventNumberAlg_check_install_runtime_makefile) RunEventNumberAlg_check_install_runtime
else
$(cmt_local_RunEventNumberAlg_check_install_runtime_makefile) : $(RunEventNumberAlg_check_install_runtime_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_RunEventNumberAlg_check_install_runtime) ] || \
	  [ ! -f $(cmt_final_setup_RunEventNumberAlg_check_install_runtime) ] || \
	  $(not_RunEventNumberAlg_check_install_runtime_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building RunEventNumberAlg_check_install_runtime.make"; \
	  $(cmtexe) -tag=$(tags) $(RunEventNumberAlg_check_install_runtime_extratags) build constituent_config -out=$(cmt_local_RunEventNumberAlg_check_install_runtime_makefile) RunEventNumberAlg_check_install_runtime; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_RunEventNumberAlg_check_install_runtime_makefile) : $(RunEventNumberAlg_check_install_runtime_dependencies) build_library_links
	$(echo) "(constituents.make) Building RunEventNumberAlg_check_install_runtime.make"; \
	  $(cmtexe) -f=$(bin)RunEventNumberAlg_check_install_runtime.in -tag=$(tags) $(RunEventNumberAlg_check_install_runtime_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_RunEventNumberAlg_check_install_runtime_makefile) RunEventNumberAlg_check_install_runtime
else
$(cmt_local_RunEventNumberAlg_check_install_runtime_makefile) : $(RunEventNumberAlg_check_install_runtime_dependencies) $(cmt_build_library_linksstamp) $(bin)RunEventNumberAlg_check_install_runtime.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_RunEventNumberAlg_check_install_runtime) ] || \
	  [ ! -f $(cmt_final_setup_RunEventNumberAlg_check_install_runtime) ] || \
	  $(not_RunEventNumberAlg_check_install_runtime_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building RunEventNumberAlg_check_install_runtime.make"; \
	  $(cmtexe) -f=$(bin)RunEventNumberAlg_check_install_runtime.in -tag=$(tags) $(RunEventNumberAlg_check_install_runtime_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_RunEventNumberAlg_check_install_runtime_makefile) RunEventNumberAlg_check_install_runtime; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(RunEventNumberAlg_check_install_runtime_extratags) build constituent_makefile -out=$(cmt_local_RunEventNumberAlg_check_install_runtime_makefile) RunEventNumberAlg_check_install_runtime

RunEventNumberAlg_check_install_runtime :: $(RunEventNumberAlg_check_install_runtime_dependencies) $(cmt_local_RunEventNumberAlg_check_install_runtime_makefile) dirs RunEventNumberAlg_check_install_runtimedirs
	$(echo) "(constituents.make) Starting RunEventNumberAlg_check_install_runtime"
	@if test -f $(cmt_local_RunEventNumberAlg_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_RunEventNumberAlg_check_install_runtime_makefile) RunEventNumberAlg_check_install_runtime; \
	  fi
#	@$(MAKE) -f $(cmt_local_RunEventNumberAlg_check_install_runtime_makefile) RunEventNumberAlg_check_install_runtime
	$(echo) "(constituents.make) RunEventNumberAlg_check_install_runtime done"

clean :: RunEventNumberAlg_check_install_runtimeclean

RunEventNumberAlg_check_install_runtimeclean :: $(RunEventNumberAlg_check_install_runtimeclean_dependencies) ##$(cmt_local_RunEventNumberAlg_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting RunEventNumberAlg_check_install_runtimeclean"
	@-if test -f $(cmt_local_RunEventNumberAlg_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_RunEventNumberAlg_check_install_runtime_makefile) RunEventNumberAlg_check_install_runtimeclean; \
	fi
	$(echo) "(constituents.make) RunEventNumberAlg_check_install_runtimeclean done"
#	@-$(MAKE) -f $(cmt_local_RunEventNumberAlg_check_install_runtime_makefile) RunEventNumberAlg_check_install_runtimeclean

##	  /bin/rm -f $(cmt_local_RunEventNumberAlg_check_install_runtime_makefile) $(bin)RunEventNumberAlg_check_install_runtime_dependencies.make

install :: RunEventNumberAlg_check_install_runtimeinstall

RunEventNumberAlg_check_install_runtimeinstall :: $(RunEventNumberAlg_check_install_runtime_dependencies) $(cmt_local_RunEventNumberAlg_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting install RunEventNumberAlg_check_install_runtime"
	@-$(MAKE) -f $(cmt_local_RunEventNumberAlg_check_install_runtime_makefile) install
	$(echo) "(constituents.make) install RunEventNumberAlg_check_install_runtime done"

uninstall : RunEventNumberAlg_check_install_runtimeuninstall

$(foreach d,$(RunEventNumberAlg_check_install_runtime_dependencies),$(eval $(d)uninstall_dependencies += RunEventNumberAlg_check_install_runtimeuninstall))

RunEventNumberAlg_check_install_runtimeuninstall : $(RunEventNumberAlg_check_install_runtimeuninstall_dependencies) ##$(cmt_local_RunEventNumberAlg_check_install_runtime_makefile)
	$(echo) "(constituents.make) Starting uninstall RunEventNumberAlg_check_install_runtime"
	@if test -f $(cmt_local_RunEventNumberAlg_check_install_runtime_makefile); then \
	  $(MAKE) -f $(cmt_local_RunEventNumberAlg_check_install_runtime_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_RunEventNumberAlg_check_install_runtime_makefile) uninstall
	$(echo) "(constituents.make) uninstall RunEventNumberAlg_check_install_runtime done"

remove_library_links :: RunEventNumberAlg_check_install_runtimeuninstall

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ RunEventNumberAlg_check_install_runtime"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ RunEventNumberAlg_check_install_runtime done"
endif

#-- end of constituent ------

clean :: remove_library_links

remove_library_links ::
	@echo "Removing library links"; \
	  $(remove_library_links); \

makefilesclean ::
	@/bin/rm -f checkuses

#	/bin/rm -f *.make*

clean :: makefilesclean

binclean :: clean
	if test ! "$(bin)" = "./"; then /bin/rm -rf $(bin); fi

