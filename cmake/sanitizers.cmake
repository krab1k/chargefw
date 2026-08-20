function(chargefw_enable_sanitizer target sanitizer)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(${target}
      PUBLIC
        -fsanitize=${sanitizer}
        -fno-omit-frame-pointer
    )

    target_link_options(${target}
      PUBLIC
        -fsanitize=${sanitizer}
    )
  endif()
endfunction()
