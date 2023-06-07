process CHECK_MAX_CHANNEL_SIZE {
    executor 'local'

    input:
        val(channel_size) //Output of Channel.count()
        val(max_size) //Integer for the max size allowed
        val(channel_name) //Channel name to be displayed in error log

    exec:
    if (channel_size > max_size) {
        println "ERROR: Channel size for $channel_name is $channel_size, but maximum allowed is $max_size"
        exit 1
    }
}

process CHECK_CHANNEL_SIZE {
    executor 'local'

    input:
        val(channel_size) //Output of Channel.count()
        val(expected_size) //Integer for the expected size
        val(channel_name) //Channel name to be displayed in error log

    exec:
    if (channel_size != expected_size) {
        println "ERROR: Channel size for $channel_name is $channel_size, but $expected_size was expected"
        exit 1
    }
}