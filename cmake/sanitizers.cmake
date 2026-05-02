function(chargefw_enable_sanitizers target)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(${target}
      PUBLIC
        -fsanitize=address,undefined
        -fno-omit-frame-pointer
    )

    target_link_options(${target}
      PUBLIC
        -fsanitize=address,undefined
    )
  endif()
endfunction()
